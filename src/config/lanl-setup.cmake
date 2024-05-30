#-----------------------------*-cmake-*----------------------------------------#
# file   src/config/lanl-setup.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Tuesday, Aug 14, 2018, 15:24 pm
# brief  Setup targeting LANL HPC
# note   Copyright (C) 2018 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

if( ${SITENAME} MATCHES "sn" OR ${SITENAME} MATCHES "ic" OR ${SITENAME} MATCHES "fi")
   set( SITENAME "Snow")
   set(SNOW_NODE TRUE)
elseif( ${SITENAME} MATCHES "ba-")
  set( SITENAME "Badger")
  set(BADGER_NODE TRUE)
elseif( ${SITENAME} MATCHES "tt-")
  set( SITENAME "Trinitite")
  set(TRINITITE_NODE TRUE)
elseif( ${SITENAME} MATCHES "tr-")
  set( SITENAME "Trinity")
  set(TRINITY_NODE TRUE)
elseif( ${SITENAME} MATCHES "ch-")
  set( SITENAME "Chicoma")
  set(CHICOMA_NODE TRUE)
elseif( ${SITENAME} MATCHES "ve-")
  set( SITENAME "Venado")
  set(VENADO_NODE TRUE)
elseif( ${SITENAME} MATCHES "ro-")
  set( SITENAME "Rocinante")
  set(ROCINANTE_NODE TRUE)
elseif( ${SITENAME} MATCHES "xr-")
  set( SITENAME "Crossroads")
  set(CROSSROADS_NODE TRUE)
elseif( ${SITENAME} MATCHES "ve-")
  set( SITENAME "Venado")
  set(VENADO_NODE TRUE)
elseif( ${SITENAME} MATCHES "tr-")
  set( SITENAME "Trinity")
  set(TRINITY_NODE TRUE)
elseif( ${SITENAME} MATCHES "ccscs[0-9]+" )
  # do nothing (keep the fullname)
  set(TRINITY_NODE TRUE)
elseif( ${SITENAME} MATCHES "nid[0-9]+" )
  set( SITENAME "CRAY_BACKEND")
  set(CRAY_BACKEND_NODE TRUE)
endif()

#------------------------------------------------------------------------------#
function(report_lanl_hpc_features)

  if( SNOW_NODE )
    message( STATUS "LANL option     : '-DSNOW_NODE'" )
  elseif( BADGER_NODE )
    message( STATUS "LANL option     : '-DBADGER_NODE'" )
  elseif( TRINITITE_NODE )
    message( STATUS "LANL option     : '-DTRINITITE_NODE'" )
  elseif( TRINITY_NODE )
    message( STATUS "LANL option     : '-DTRINITY_NODE'" )
  endif()

endfunction()

#------------------------------------------------------------------------------#
# End lanl-setup.cmake
#------------------------------------------------------------------------------#
