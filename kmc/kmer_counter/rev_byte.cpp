#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.1.1
  Date   : 2015-01-22
*/

#include "rev_byte.h"

uchar CRev_byte::lut[256];
CRev_byte::_si CRev_byte::_init;