from setuptools import setup, find_packages

setup(
    name='cfrna-source-tracing',
    version='5.2.0',
    description='cfRNA source tracing toolkit with paper-grade benchmark export',
    packages=find_packages(include=['app', 'app.*', 'benchmark', 'benchmark.*', 'core', 'core.*', 'data', 'data.*', 'reporting', 'reporting.*']),
    py_modules=['benchmark_runner', 'data_processor', 'database_init', 'signature_builder', 'source_tracing', 'source_tracing_v2', 'cfrna_tracing_app', 'cli'],
    include_package_data=True,
    install_requires=[
        'streamlit>=1.28.0',
        'pandas>=2.0.0',
        'numpy>=1.24.0',
        'plotly>=5.17.0',
        'scikit-learn>=1.3.0',
        'scipy>=1.11.0',
        'openpyxl>=3.1.0',
        'matplotlib>=3.8.0',
    ],
    entry_points={
        'console_scripts': [
            'cfrna-tracing=cli:main',
        ]
    },
)
