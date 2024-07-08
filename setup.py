from setuptools import setup


long_description = open('README.md', encoding='utf-8').read()

setup(name='LorBin',
      version='0.1.0',
      python_requires='>=3.10',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      #url='https://github.com/'
      #license='MIT',
      packages = ['lorbin','lorbin/model'],
      package_data={
          'lorbin': ['*.hmm', 'model/*.pt','model/*.csv']},
      zip_safe=False,
      entry_points={
          'console_scripts': ['LorBin=lorbin.lorbin:main',
                                ],
      }
)

