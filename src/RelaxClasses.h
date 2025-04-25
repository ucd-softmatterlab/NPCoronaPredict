#ifndef REL_CLASS__H__
#define REL_CLASS__H__

#include "StringFormat.h"

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <unordered_map>

class bond  {
public:
    int i;
    int j;
    double length;
    double bondk;

    bond( int iIn, int jIn, double lengthIn, double kIn){
     i = iIn;
     j = jIn;
     length = lengthIn;
     bondk = kIn; 
    }

};


#endif
