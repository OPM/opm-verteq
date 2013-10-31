#
# spec file for package opm-verteq
#

%define tag rc2

Name:           opm-verteq
Version:        2013.10
Release:        0
Summary:        Open Porous Media - vertical equilibrium library
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires:  blas-devel lapack-devel dune-common-devel
BuildRequires:  git suitesparse-devel cmake28 doxygen bc
%{?el5:BuildRequires: gcc44 gcc44-c++}
%{!?el5:BuildRequires: gcc gcc-c++}
%{?el5:BuildRequires: boost141-devel}
%{!?el5:BuildRequires: boost-devel}
BuildRequires:  tinyxml-devel dune-istl-devel opm-core-devel ert.ecl-devel
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package -n libopm-verteq0
Summary:        Open Porous Media - vertical equilibrium library
Group:          System/Libraries
%{?el5:BuildArch: %{_arch}}

%description -n libopm-verteq0
The Open Porous Media (OPM) initiative provides a set of open-source tools centered around the simulation of flow and transport of fluids in porous media. The goal of the initiative is to establish a sustainable environment for the development of an efficient and well-maintained software suite.

%package devel
Summary:        Development and header files for opm-verteq
Group:          Development/Libraries/C and C++
%{?el5:BuildArch: %{_arch}}

%description devel
This package contains the development and header files for opm-verteq

%package doc
Summary:        Documentation files for opm-verteq
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-verteq

%{?el5:
%package debuginfo
Summary:        Debug info in opm-verteq
Group:          Scientific
Requires:       libopm-verteq0 = %{version}
BuildArch: 	%{_arch}

%description debuginfo
This package contains the debug symbols for opm-verteq
}

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF %{?el5:-DCMAKE_CXX_COMPILER=g++44 -DCMAKE_C_COMPILER=gcc44 -DBOOST_LIBRARYDIR=%{_libdir}/boost141 -DBOOST_INCLUDEDIR=/usr/include/boost141}
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%files doc
%{_docdir}/*

%files -n libopm-verteq0
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*

%{?el5:
%files debuginfo
/usr/lib/debug/%{_libdir}/*.so*.debug
}
