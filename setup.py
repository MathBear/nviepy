from distutils.core import setup

setup(name='nviepy',
      version='0.1',
      package_dir={'nviepy': ''},
      packages=['nviepy'],
      #package_data={'' : ['examples/*.py']},
      author=['Vincenzo Schiano'],
      author_email=['vin.schianodicola@studenti.unina.it'],
      url='https://github.com/MathBear/nviepy',
      description='Numerical VIE order conditions creator',
      license='modified BSD',
      requires=['numpy','sympy'],
      )
