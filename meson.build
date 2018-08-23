project('coastal-engineering')

# Importing modules
python = import('python3')

# Make sure Python is installed and found
if not python.find_python().found()
    error('No valid python3 binary found')
endif

PYTHON_DIR = join_paths(get_option('prefix'), python.sysconfig_path('purelib'))

# Install the Python package
install_subdir(
    'coastal_engineering',
    install_dir: PYTHON_DIR
)

# Install the executable file
bin_config = configuration_data()
bin_config.set('pythondir', PYTHON_DIR)
configure_file(
    input: 'coastal-engineering.in',
    output: 'coastal-engineering',
    configuration: bin_config,
    install_dir: get_option('bindir')
)