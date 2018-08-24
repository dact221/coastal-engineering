project('coastal-engineering')

# Importing modules
python = import('python3')

# Make sure Python is installed and found
if not python.find_python().found()
    error('No valid python3 binary found')
endif

# Check Python version
if python.language_version().version_compare('<3.6')
    error('You need Python 3.6')
endif

PYTHON_DIR = join_paths(get_option('prefix'), python.sysconfig_path('purelib'))

# Install the Python package
install_subdir(
    'coastal_engineering',
    exclude_directories: '__pycache__',
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

# Install the builddir executable
local_config = configuration_data()
local_config.set('pythondir', meson.source_root())
configure_file(
    input: 'coastal-engineering.in',
    output: 'local-coastal-engineering',
    configuration: local_config
)
