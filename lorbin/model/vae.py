import torch
from torch import nn
from typing import Optional, IO, Union
from pathlib import Path
from torch import Tensor
from torch.nn.functional import softmax
from math import log
from torch.utils.data.dataset import TensorDataset
from torch.utils.data import DataLoader
from torch.optim import Adam
import numpy as np
import datetime
import logging

def zscore(array: np.ndarray, axis: Optional[int] = None) -> np.ndarray:
    """Calculates zscore for an array. A cheap copy of scipy.stats.zscore.

    Inputs:
        array: Numpy array to be normalized
        axis: Axis to operate across [None = entrie array]
        inplace: Do not create new array, change input array [False]

    Output:
        If inplace is True: Input numpy array
        else: New normalized Numpy-array
    """

    if axis is not None and (axis >= array.ndim or axis < 0):
        raise np.AxisError(str(axis))


    mean = array.mean(axis=axis)
    std = array.std(axis=axis)

    if axis is None:
        if std == 0:
            std = 1  # prevent divide by zero

    else:
        std[std == 0.0] = 1  # prevent divide by zero
        shape = tuple(dim if ax != axis else 1 for ax, dim in enumerate(array.shape))
        mean.shape, std.shape = shape, shape

    return (array - mean) / std




def normalize(rpkm: np.ndarray,tnf: np.ndarray,lengths: np.ndarray):
    if not isinstance(rpkm, np.ndarray) or not isinstance(tnf, np.ndarray):
        raise ValueError("TNF and RPKM must be Numpy arrays")

    if len(rpkm) != len(tnf) or len(tnf) != len(lengths):
        raise ValueError("Lengths of RPKM, TNF and lengths arrays must be the same")

    if not (rpkm.dtype == tnf.dtype == np.float32):
        raise ValueError("TNF and RPKM must be Numpy arrays of dtype float32")

    rpkm = rpkm.copy()
    tnf = tnf.copy()

    sample_depths_sum = rpkm.sum(axis=0)
    if np.any(sample_depths_sum == 0):
        raise ValueError(
            "One or more samples have zero depth in all sequences, so cannot be depth normalized"
        )
    rpkm *= 1_000_000 / sample_depths_sum

    zero_tnf = tnf.sum(axis=1) == 0
    smallest_index = np.argmax(zero_tnf)
    if zero_tnf[smallest_index]:
        raise ValueError(
            f"TNF row at index {smallest_index} is all zeros. "
            + "This implies that the sequence contained no 4-mers of A, C, G, T or U, "
            + "making this sequence uninformative. This is probably a mistake. "
            + "Verify that the sequence contains usable information (e.g. is not all N's)"
        )

    total_abundance = rpkm.sum(axis=1)

    # Normalize rpkm to sum to 1
    n_samples = rpkm.shape[1]
    zero_total_abundance = total_abundance == 0
    rpkm[zero_total_abundance] = 1 / n_samples
    nonzero_total_abundance = total_abundance.copy()
    nonzero_total_abundance[zero_total_abundance] = 1.0
    rpkm /= nonzero_total_abundance.reshape((-1, 1))

    # Normalize TNF and total abundance to make SSE loss work better
    total_abundance = np.log(total_abundance.clip(min=0.001))
    total_abundance=zscore(total_abundance)
    tnf=zscore(tnf, axis=0)
    total_abundance.shape = (len(total_abundance), 1)

    # Create weights
    lengths = (lengths).astype(np.float32)
    weights = np.log(lengths).astype(np.float32) - 5.0
    weights[weights < 2.0] = 2.0
    weights *= len(weights) / weights.sum()
    weights.shape = (len(weights), 1)

    return rpkm, tnf, total_abundance, weights

def make_dataloader(rpkm:np.ndarray,tnf:np.ndarray,total_abundance:np.ndarray,
                    weights:np.ndarray,batchsize:int=256)->DataLoader:
    depthstensor = torch.from_numpy(rpkm)  # this is a no-copy operation
    tnftensor = torch.from_numpy(tnf)
    total_abundance_tensor = torch.from_numpy(total_abundance)
    weightstensor = torch.from_numpy(weights)

    dataset = TensorDataset(depthstensor, tnftensor, total_abundance_tensor, weightstensor)
    dataloader = DataLoader(dataset=dataset,batch_size=batchsize,drop_last=True,shuffle=True)
    return dataloader



class VAE(nn.Module):
    """Variational autoencoder, subclass of torch.nn.Module.

    Instantiate with:
        nsamples: Number of samples in abundance matrix
        nhiddens: list of n_neurons in the hidden layers [None=Auto]
        nlatent: Number of neurons in the latent layer [32]
        alpha: Approximate starting TNF/(CE+TNF) ratio in loss. [None = Auto]
        beta: Multiply KLD by the inverse of this value [200]
        dropout: Probability of dropout on forward pass [0.2]
        cuda: Use CUDA (GPU accelerated training) [False]

    vae.trainmodel(dataloader, nepochs batchsteps, lrate, logfile, modelfile)
        Trains the model, returning None

    vae.encode(self, data_loader):
        Encodes the data in the data loader and returns the encoded matrix.

    If alpha or dropout is None and there is only one sample, they are set to
    0.99 and 0.0, respectively
    """

    def __init__(
        self,
        nsamples: int,
        nhiddens: Optional[list[int]] = None,
        nlatent: int = 32,
        alpha: Optional[float] = None,
        beta: float = 200.0,
        dropout: Optional[float] = 0.2,
        cuda: bool = False,
        seed: int = 0,
    ):
        if nlatent < 1:
            raise ValueError(f"Minimum 1 latent neuron, not {nlatent}")

        if nsamples < 1:
            raise ValueError(f"nsamples must be > 0, not {nsamples}")

        # If only 1 sample, we weigh alpha and nhiddens differently
        if alpha is None:
            alpha = 0.15 if nsamples > 1 else 0.50

        if nhiddens is None:
            nhiddens = [512, 512] if nsamples > 1 else [256, 256]

        if dropout is None:
            dropout = 0.2 if nsamples > 1 else 0.0

        if any(i < 1 for i in nhiddens):
            raise ValueError(f"Minimum 1 neuron per layer, not {min(nhiddens)}")

        if beta <= 0:
            raise ValueError(f"beta must be > 0, not {beta}")

        if not (0 < alpha < 1):
            raise ValueError(f"alpha must be 0 < alpha < 1, not {alpha}")

        if not (0 <= dropout < 1):
            raise ValueError(f"dropout must be 0 <= dropout < 1, not {dropout}")

        torch.manual_seed(seed)
        super(VAE, self).__init__()

        # Initialize simple attributes
        self.usecuda = cuda
        self.nsamples = nsamples
        self.ntnf = 136
        self.alpha = alpha
        self.beta = beta
        self.nhiddens = nhiddens
        self.nlatent = nlatent
        self.dropout = dropout

        # Initialize lists for holding hidden layers
        self.encoderlayers = nn.ModuleList()
        self.encodernorms = nn.ModuleList()
        self.decoderlayers = nn.ModuleList()
        self.decodernorms = nn.ModuleList()

        # Add all other hidden layers
        for nin, nout in zip(
            # + 1 for the total abundance
            [self.nsamples + self.ntnf + 1] + self.nhiddens,
            self.nhiddens,
        ):
            self.encoderlayers.append(nn.Linear(nin, nout))
            self.encodernorms.append(nn.BatchNorm1d(nout))

        # Latent layers
        self.mu = nn.Linear(self.nhiddens[-1], self.nlatent)

        # Add first decoding layer
        for nin, nout in zip([self.nlatent] + self.nhiddens[::-1], self.nhiddens[::-1]):
            self.decoderlayers.append(nn.Linear(nin, nout))
            self.decodernorms.append(nn.BatchNorm1d(nout))

        # Reconstruction (output) layer. + 1 for the total abundance
        self.outputlayer = nn.Linear(self.nhiddens[0], self.nsamples + self.ntnf + 1)

        # Activation functions
        self.relu = nn.LeakyReLU()
        self.softplus = nn.Softplus()
        self.dropoutlayer = nn.Dropout(p=self.dropout)

        if cuda:
            self.cuda()

    def encode(self, tensor: Tensor) -> Tensor:
        tensors = list()

        # Hidden layers
        for encoderlayer, encodernorm in zip(self.encoderlayers, self.encodernorms):
            tensor = encodernorm(self.dropoutlayer(self.relu(encoderlayer(tensor))))
            tensors.append(tensor)

        # Latent layers
        mu = self.mu(tensor)

        # Note: We ought to also compute logsigma here, but we had a bug in the original
        # implementation of Vamb where logsigma was fixed to zero, so we just remove it.

        return mu

    # sample with gaussian noise
    def reparameterize(self, mu: Tensor) -> Tensor:
        epsilon = torch.randn(mu.size(0), mu.size(1))
        if self.usecuda:
            epsilon = epsilon.cuda()
        epsilon.requires_grad = True
        latent = mu + epsilon
        return latent

    def decode(self, tensor: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        tensors = list()

        for decoderlayer, decodernorm in zip(self.decoderlayers, self.decodernorms):
            tensor = decodernorm(self.dropoutlayer(self.relu(decoderlayer(tensor))))
            tensors.append(tensor)

        reconstruction = self.outputlayer(tensor)

        # Decompose reconstruction to depths and tnf signal
        depths_out = reconstruction.narrow(1, 0, self.nsamples)
        tnf_out = reconstruction.narrow(1, self.nsamples, self.ntnf)
        abundance_out = reconstruction.narrow(1, self.nsamples + self.ntnf, 1)

        depths_out = softmax(depths_out, dim=1)

        return depths_out, tnf_out, abundance_out

    def forward(
        self, depths: Tensor, tnf: Tensor, abundance: Tensor
    ) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        tensor = torch.cat((depths, tnf, abundance), 1)
        mu = self.encode(tensor)
        latent = self.reparameterize(mu)
        depths_out, tnf_out, abundance_out = self.decode(latent)

        return depths_out, tnf_out, abundance_out, mu

    def calc_loss(
        self,
        depths_in: Tensor,
        depths_out: Tensor,
        tnf_in: Tensor,
        tnf_out: Tensor,
        abundance_in: Tensor,
        abundance_out: Tensor,
        mu: Tensor,
        weights: Tensor,
    ) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        ab_sse = (abundance_out - abundance_in).pow(2).sum(dim=1)
        # Add 1e-9 to depths_out to avoid numerical instability.
        ce = -((depths_out + 1e-9).log() * depths_in).sum(dim=1)
        sse = (tnf_out - tnf_in).pow(2).sum(dim=1)
        kld = 0.5 * (mu.pow(2)).sum(dim=1)

        # Avoid having the denominator be zero
        if self.nsamples == 1:
            ce_weight = 0.0
        else:
            ce_weight = ((1 - self.alpha) * (self.nsamples - 1)) / (
                self.nsamples * log(self.nsamples)
            )

        ab_sse_weight = (1 - self.alpha) * (1 / self.nsamples)
        sse_weight = self.alpha / self.ntnf
        kld_weight = 1 / (self.nlatent * self.beta)
        weighed_ab = ab_sse * ab_sse_weight
        weighed_ce = ce * ce_weight
        weighed_sse = sse * sse_weight
        weighed_kld = kld * kld_weight
        reconstruction_loss = weighed_ce + weighed_ab + weighed_sse
        loss = (reconstruction_loss + weighed_kld) * weights

        return (
            loss.mean(),
            weighed_ab.mean(),
            weighed_ce.mean(),
            weighed_sse.mean(),
            weighed_kld.mean(),
        )

    def trainepoch(
        self,
        logger,
        dataloader: DataLoader,
        epoch: int,
        optimizer,
        batchsteps: list[int],
    ) -> DataLoader[tuple[Tensor, Tensor, Tensor]]:
        self.train()

        epoch_loss = 0.0
        epoch_kldloss = 0.0
        epoch_sseloss = 0.0
        epoch_celoss = 0.0
        epoch_absseloss = 0.0

        if epoch in batchsteps:
            dataloader = DataLoader(dataloader.dataset,batch_size=dataloader.batch_size*2)

        for depths_in, tnf_in, abundance_in, weights in dataloader:
            depths_in.requires_grad = True
            tnf_in.requires_grad = True
            abundance_in.requires_grad = True

            if self.usecuda:
                depths_in = depths_in.cuda()
                tnf_in = tnf_in.cuda()
                abundance_in = abundance_in.cuda()
                weights = weights.cuda()

            optimizer.zero_grad()

            depths_out, tnf_out, abundance_out, mu = self(
                depths_in, tnf_in, abundance_in
            )

            loss, ab_sse, ce, sse, kld = self.calc_loss(
                depths_in,
                depths_out,
                tnf_in,
                tnf_out,
                abundance_in,
                abundance_out,
                mu,
                weights,
            )

            loss.backward()
            optimizer.step()

            epoch_loss += loss.data.item()
            epoch_kldloss += kld.data.item()
            epoch_sseloss += sse.data.item()
            epoch_celoss += ce.data.item()
            epoch_absseloss += ab_sse.data.item()


        
        logger.info("\tTime: {}\tEpoch: {:>3}  Loss: {:.5e}  CE: {:.5e}  AB: {:.5e}  SSE: {:.5e}  KLD: {:.5e}  Batchsize: {}".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                epoch + 1,
                epoch_loss / len(dataloader),
                epoch_celoss / len(dataloader),
                epoch_absseloss / len(dataloader),
                epoch_sseloss / len(dataloader),
                epoch_kldloss / len(dataloader),
                dataloader.batch_size,
            ))

        self.eval()
        return dataloader

    def get_latent(self, data_loader) -> np.ndarray:
        """Encode a data loader to a latent representation with VAE

        Input: data_loader: As generated by train_vae

        Output: A (n_contigs x n_latent) Numpy array of latent repr.
        """

        self.eval()


        depths_array, _, _, _ = data_loader.dataset.tensors
        length = len(depths_array)

        # We make a Numpy array instead of a Torch array because, if we create
        # a Torch array, then convert it to Numpy, Numpy will believe it doesn't
        # own the memory block, and array resizes will not be permitted.
        latent = np.empty((length, self.nlatent), dtype=np.float32)

        row = 0
        with torch.no_grad():
            for depths, tnf, ab, _ in data_loader:
                # Move input to GPU if requested
                if self.usecuda:
                    depths = depths.cuda()
                    tnf = tnf.cuda()
                    ab = ab.cuda()

                # Evaluate
                _, _, _, mu = self(depths, tnf, ab)

                if self.usecuda:
                    mu = mu.cpu()

                latent[row : row + len(mu)] = mu
                row += len(mu)
        print(row,length)
        assert row == length
        return latent

    def save(self, filehandle):
        """Saves the VAE to a path or binary opened file. Load with VAE.load

        Input: Path or binary opened filehandle
        Output: None
        """
        state = {
            "nsamples": self.nsamples,
            "alpha": self.alpha,
            "beta": self.beta,
            "dropout": self.dropout,
            "nhiddens": self.nhiddens,
            "nlatent": self.nlatent,
            "state": self.state_dict(),
        }

        torch.save(state, filehandle)

    @classmethod
    def load(
        cls, path: Union[IO[bytes], str], cuda: bool = False, evaluate: bool = True
    ):
        """Instantiates a VAE from a model file.

        Inputs:
            path: Path to model file as created by functions VAE.save or
                  VAE.trainmodel.
            cuda: If network should work on GPU [False]
            evaluate: Return network in evaluation mode [True]

        Output: VAE with weights and parameters matching the saved network.
        """

        # Forcably load to CPU even if model was saves as GPU model
        dictionary = torch.load(path, map_location=lambda storage, loc: storage)

        nsamples = dictionary["nsamples"]
        alpha = dictionary["alpha"]
        beta = dictionary["beta"]
        dropout = dictionary["dropout"]
        nhiddens = dictionary["nhiddens"]
        nlatent = dictionary["nlatent"]
        state = dictionary["state"]

        vae = cls(nsamples, nhiddens, nlatent, alpha, beta, dropout, cuda)
        vae.load_state_dict(state)

        if cuda:
            vae.cuda()

        if evaluate:
            vae.eval()

        return vae

    def trainmodel(
        self,
        logger,
        outdir,
        dataloader: DataLoader[tuple[Tensor, Tensor, Tensor]],
        nepochs: int = 500,
        lrate: float = 1e-3,
        batchsteps: Optional[list[int]] = [25, 75, 150, 300]
    ):
        """Train the autoencoder from depths array and tnf array.

        Inputs:
            logger: 运行日记
            dataloader: DataLoader made by make_dataloader
            nepochs: Train for this many epochs before encoding [500]
            lrate: Starting learning rate for the optimizer [0.001]
            batchsteps: None or double batchsize at these epochs [25, 75, 150, 300]
            logfile: Print status updates to this file if not None [None]

        Output: None
        """

        if lrate < 0:
            raise ValueError(f"Learning rate must be positive, not {lrate}")

        if nepochs < 1:
            raise ValueError("Minimum 1 epoch, not {nepochs}")

        if batchsteps is None:
            batchsteps_set: set[int] = set()
        else:
            # First collect to list in order to allow all element types, then check that
            # they are integers
            batchsteps = list(batchsteps)
            if not all(isinstance(i, int) for i in batchsteps):
                raise ValueError("All elements of batchsteps must be integers")
            if max(batchsteps, default=0) >= nepochs:
                raise ValueError("Max batchsteps must not equal or exceed nepochs")
            last_batchsize = dataloader.batch_size * 2 ** len(batchsteps)
            if len(dataloader.dataset) < last_batchsize:  # type: ignore
                raise ValueError(
                    f"Last batch size of {last_batchsize} exceeds dataset length "
                    f"of {len(dataloader.dataset)}. "  # type: ignore
                    "This means you have too few contigs left after filtering to train. "
                    "It is not adviced to run Vamb with fewer than 10,000 sequences "
                    "after filtering. "
                    "Please check the Vamb log file to see where the sequences were "
                    "filtered away, and verify BAM files has sensible content."
                )
            batchsteps_set = set(batchsteps)

        # Get number of features
        # Following line is un-inferrable due to typing problems with DataLoader
        ncontigs, nsamples = dataloader.dataset.tensors[0].shape  # type: ignore
        optimizer = Adam(self.parameters(), lr=lrate)

        logger.info("\tNetwork properties:")
        logger.info(f"\tCUDA: {self.usecuda}")
        logger.info(f"\tAlpha:{self.alpha}")
        logger.info(f"\tBeta:{self.beta}")
        logger.info(f"\tDropout:{self.dropout}")
        logger.info(f"\tN hidden:f{ ', '.join(map(str, self.nhiddens))}")
        logger.info(f"\tN latent:{self.nlatent}")
        logger.info("\n\tTraining properties:")
        logger.info(f"\tN epochs:{nepochs}")
        logger.info(f"\tStarting batch size:{dataloader.batch_size}")
        batchsteps_string = (
            ", ".join(map(str, sorted(batchsteps_set)))
            if batchsteps_set
             else "None"
        )
        logger.info(f"\tBatchsteps:{batchsteps_string}")
        logger.info(f"\tLearning rate:{lrate}")
        logger.info(f"\tN sequences: {ncontigs}")
        logger.info(f"\tN samples:{self.nsamples}")

        # Train
        for epoch in range(nepochs):
            dataloader = self.trainepoch(logger, dataloader, epoch, optimizer, sorted(batchsteps_set))

        # Save weights - Lord forgive me, for I have sinned when catching all exceptions
        try:
            self.save(f"{outdir}/model.pt")
        except:
            pass

        return None
