from setuptools import setup

setup(name='limathpy',
      version='0.2',
      description='A Python package for undergraduate math',
      url='http://github.com/YessiRocha/pylima',
      author='Yessica Rocha',
      author_email='yesi.bundo@gmail.com',
      license='MIT',
      packages=['limathpy'],
      install_requires=['sympy', 'numpy', 'matplotlib'],
      zip_safe=False)
