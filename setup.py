from setuptools import setup
pname='halosky'
setup(name=pname,
      version='0.1',
      description='Random SZ maps',
      url='http://github.com/marcelo-alvarez/halosky',
      author='Marcelo Alvarez',
      license='MIT',
      packages=['halosky'],
      package_data={
        pname: ["data/*"]
      },
      zip_safe=False)
