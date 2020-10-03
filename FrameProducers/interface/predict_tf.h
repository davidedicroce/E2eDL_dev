#ifndef predict_tf_h
#define predict_tf_h

#include <vector>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "E2eDL/DataFormats/interface/FrameCollections.h"
#include <iostream>
#include <fstream>
using namespace std;

e2e::Frame2D predict_tf(e2e::Frame4D&, string, string, string);

#endif
