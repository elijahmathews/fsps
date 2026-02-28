doxygen = find_program('doxygen', required : true)
sphinx_build = find_program('sphinx-build', required : true)
pip = find_program('pip', required : true)
sh = find_program('sh', required : true)

doxygen_output_dir = join_paths(meson.current_build_dir(), '_build', 'doxygen')
doxygen_xml_dir = join_paths(doxygen_output_dir, 'xml')
sphinx_html_dir = join_paths(meson.current_build_dir(), '_build', 'html')

doxy_cfg = configuration_data()
doxy_cfg.set('PROJECT_VERSION', meson.project_version())
doxy_cfg.set('DOXYGEN_OUTPUT_DIR', doxygen_output_dir)
doxy_cfg.set('FSPS_SRC_DIR', join_paths(meson.project_source_root(), 'src'))
doxy_cfg.set('FSPS_INCLUDE_DIR', join_paths(meson.project_source_root(), 'include'))

configured_doxyfile = configure_file(
    input : 'Doxyfile.in',
    output : 'Doxyfile',
    configuration : doxy_cfg,
)

doxygen_xml_target = custom_target(
    'docs_doxygen_xml',
    input : configured_doxyfile,
    output : 'doxygen-xml.stamp',
    command : [
        sh,
        '-c',
        'set -eu; mkdir -p "$2"; "$1" "$3"; touch "$4"',
        'docs-doxygen',
        doxygen,
        doxygen_output_dir,
        '@INPUT@',
        '@OUTPUT@',
    ],
    build_by_default : true,
    console : true,
)

sphinx_html_target = custom_target(
    'docs_sphinx_html',
    input : [
        files('index.rst', 'fortran_api.rst', 'physics.rst', 'c_api.md', 'conf.py'),
        configured_doxyfile,
    ],
    output : 'docs-html.stamp',
    depends : [doxygen_xml_target],
    command : [
        sh,
        '-c',
        'set -eu; mkdir -p "$1"; DOXYGEN_XML_DIR="$2" "$3" -b html "$4" "$1"; touch "$5"',
        'docs-sphinx',
        sphinx_html_dir,
        doxygen_xml_dir,
        sphinx_build,
        meson.current_source_dir(),
        '@OUTPUT@',
    ],
    build_by_default : true,
    console : true,
)

doc_install_dir = join_paths(get_option('datadir'), 'doc', 'fsps', 'html')

meson.add_install_script(
    sh,
    '-c',
    'set -eu; src="$1"; rel="$2"; dst="${DESTDIR:-}${MESON_INSTALL_PREFIX}/$rel"; mkdir -p "$dst"; cp -a "$src"/. "$dst"/',
    'install-docs',
    sphinx_html_dir,
    doc_install_dir,
)
