from setuptools import setup

setup(name='pylima',
      version='0.1',
      description='A Python package for undergraduate math',
      url='http://github.com/YessiRocha/pylima',
      author='Yessica Rocha',
      author_email='yesi.bundo@gmail.com',
      license='MIT',
      packages=['pylima'],
      install_requires=['sympy', 'numpy', 'matplotlib'],
      zip_safe=False)
