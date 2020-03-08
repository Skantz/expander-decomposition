# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_ALL_INSTALL_TYPES "Full;Developer")
set(CPACK_BINARY_7Z "")
set(CPACK_BINARY_BUNDLE "")
set(CPACK_BINARY_CYGWIN "")
set(CPACK_BINARY_DEB "")
set(CPACK_BINARY_DRAGNDROP "")
set(CPACK_BINARY_FREEBSD "")
set(CPACK_BINARY_IFW "")
set(CPACK_BINARY_NSIS "")
set(CPACK_BINARY_NUGET "")
set(CPACK_BINARY_OSXX11 "")
set(CPACK_BINARY_PACKAGEMAKER "")
set(CPACK_BINARY_PRODUCTBUILD "")
set(CPACK_BINARY_RPM "")
set(CPACK_BINARY_STGZ "")
set(CPACK_BINARY_TBZ2 "")
set(CPACK_BINARY_TGZ "")
set(CPACK_BINARY_TXZ "")
set(CPACK_BINARY_TZ "")
set(CPACK_BINARY_WIX "")
set(CPACK_BINARY_ZIP "")
set(CPACK_BUILD_SOURCE_DIRS "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1;/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENTS_ALL "headers;library;html_documentation;bin")
set(CPACK_COMPONENTS_ALL_SET_BY_USER "TRUE")
set(CPACK_COMPONENT_BIN_DESCRIPTION "Command line utilities")
set(CPACK_COMPONENT_BIN_DISPLAY_NAME "Command line utilities")
set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "Components needed to develop software using LEMON")
set(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION "Documentation of LEMON")
set(CPACK_COMPONENT_HEADERS_DEPENDS "library")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ header files")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ headers")
set(CPACK_COMPONENT_HEADERS_GROUP "Development")
set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES "Developer;Full")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION "Doxygen generated documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")
set(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES "Full")
set(CPACK_COMPONENT_LIBRARY_DESCRIPTION "DLL and import library")
set(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME "Dynamic-link library")
set(CPACK_COMPONENT_LIBRARY_GROUP "Development")
set(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES "Developer;Full")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_GENERATOR "NSIS")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build;LEMON;ALL;/")
set(CPACK_INSTALL_PREFIX "/usr/local")
set(CPACK_MODULE_PATH "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/cmake")
set(CPACK_NSIS_CONTACT "lemon-user@lemon.cs.elte.hu")
set(CPACK_NSIS_CREATE_ICONS_EXTRA "
    CreateShortCut \"$SMPROGRAMS\\$STARTMENU_FOLDER\\Documentation.lnk\" \"$INSTDIR\\share\\doc\\index.html\"
    ")
set(CPACK_NSIS_DELETE_ICONS_EXTRA "
    !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    Delete \"$SMPROGRAMS\\$MUI_TEMP\\Documentation.lnk\"
    ")
set(CPACK_NSIS_DISPLAY_NAME "LEMON 1.3.1 LEMON")
set(CPACK_NSIS_DISPLAY_NAME_SET "TRUE")
set(CPACK_NSIS_HELP_LINK "http:\\\\lemon.cs.elte.hu")
set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\lemon.ico")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_MUI_ICON "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/cmake/nsis/lemon.ico")
set(CPACK_NSIS_MUI_UNIICON "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/cmake/nsis/uninstall.ico")
set(CPACK_NSIS_PACKAGE_NAME "LEMON 1.3.1 LEMON")
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\lemon.cs.elte.hu")
set(CPACK_OUTPUT_CONFIG_FILE "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-3.12/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LEMON - Library for Efficient Modeling and Optimization in Networks")
set(CPACK_PACKAGE_FILE_NAME "LEMON-1.3.1-Linux")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "LEMON 1.3.1")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "LEMON 1.3.1")
set(CPACK_PACKAGE_NAME "LEMON")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "EGRES")
set(CPACK_PACKAGE_VERSION "1.3.1")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/LICENSE")
set(CPACK_RESOURCE_FILE_README "/usr/share/cmake-3.12/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake-3.12/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_7Z "")
set(CPACK_SOURCE_CYGWIN "")
set(CPACK_SOURCE_GENERATOR "TBZ2;TGZ;TXZ;TZ")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/CPackSourceConfig.cmake")
set(CPACK_SOURCE_RPM "OFF")
set(CPACK_SOURCE_TBZ2 "ON")
set(CPACK_SOURCE_TGZ "ON")
set(CPACK_SOURCE_TXZ "ON")
set(CPACK_SOURCE_TZ "ON")
set(CPACK_SOURCE_ZIP "OFF")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
