from pathlib import Path
import os

project = "FSPS"
version = "3.2.0"
release = "3.2.0"

extensions = [
    "myst_parser",
    "breathe",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Markdown support via MyST for files such as doc/FSPS_C_API.md.
myst_enable_extensions = [
    "deflist",
    "fieldlist",
]

# Breathe / Doxygen bridge
_docs_dir = Path(__file__).resolve().parent
_default_xml = _docs_dir / "_build" / "doxygen" / "xml"
_doxygen_xml_dir = os.environ.get("DOXYGEN_XML_DIR", str(_default_xml))

breathe_projects = {
    "FSPS": _doxygen_xml_dir,
}
breathe_default_project = "FSPS"
breathe_domain_by_extension = {
    "f90": "fortran",
    "F90": "fortran",
    "h": "c",
    "c": "c",
}
breathe_domain_by_file_pattern = {
    "*fsps_api.f90": "fortran",
}

primary_domain = "c"
default_role = "any"

html_theme = os.environ.get("SPHINX_THEME", "furo")
