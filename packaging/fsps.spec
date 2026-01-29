Name:           fsps
Version:        3.2.0
Release:        1%{?dist}
Summary:        Flexible Stellar Population Synthesis (FSPS) command-line drivers
Requires:       libfsps%{?_isa} = %{version}-%{release}
Requires:       libfsps-data = %{version}-%{release}

License:        MIT
URL:            https://github.com/cconroy20/fsps
Source0:        %{name}-%{version}.tar.gz

BuildRequires:  gcc-gfortran
BuildRequires:  make

%description
FSPS is a stellar population synthesis code. This package provides the
command-line drivers (simple, autosps, lesssimple, spec_bin).

%package -n libfsps
Summary:        FSPS shared library and runtime data
Requires:       libgfortran
Requires:       libfsps-data = %{version}-%{release}

%description -n libfsps
Shared library for FSPS and the data files required at runtime. This is
the package typically needed by language wrappers.

%package -n libfsps-data
Summary:        FSPS runtime data files

%description -n libfsps-data
Physical libraries, templates, and tables required by FSPS at runtime.

%package -n libfsps-devel
Summary:        Development files for libfsps
Requires:       libfsps%{?_isa} = %{version}-%{release}
Provides:       pkgconfig(fsps)

%description -n libfsps-devel
Header files, unversioned linker symlink, and pkg-config metadata for libfsps.

%prep
%setup -q

%build
%make_build all shared

%install
%make_install \
    PREFIX=%{_prefix} \
    LIBDIR=%{_libdir} \
    INCLUDEDIR=%{_includedir} \
    DATADIR=%{_datadir} \
    PKGCONFIGDIR=%{_libdir}/pkgconfig \
    BINDIR=%{_bindir}

%check
FSPS_DATA_HOME=%{_builddir}/%{name}-%{version} \
%make_build check

%files
%license LICENSE
%{_bindir}/simple
%{_bindir}/lesssimple
%{_bindir}/autosps
%{_bindir}/spec_bin

%files -n libfsps
%license LICENSE
%{_libdir}/libfsps.so.*

%files -n libfsps-data
%license LICENSE
%{_datadir}/fsps/data

%files -n libfsps-devel
%license LICENSE
%{_includedir}/fsps.h
%{_libdir}/libfsps.so
%{_libdir}/pkgconfig/fsps.pc

%changelog
* Thu Jan 29 2026 FSPS Packager <packager@example.com> - 3.2.0-1
- Initial RPM spec for FSPS