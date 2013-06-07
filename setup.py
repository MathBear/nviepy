from distutils.core import setup

setup(name='nviepy',
      version='0.1',
      package_dir={'nodepy': ''},
      packages=['nodepy'],
      #package_data={'' : ['examples/*.py']},
      author=['Vincenzo Schiano'],
      author_email=['vin.schianodicola@studenti.unina.it'],
      url='http://numerics.kaust.edu.sa/nodepy/',
      description='Numerical VIE order conditions creator',
      license='BSD',
      requires=['numpy','sympy'],
      )
