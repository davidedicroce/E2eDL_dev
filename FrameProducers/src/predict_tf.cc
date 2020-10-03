//#include "ProdTutorial/ProducerTest/plugins/predict_tf.h"
#include "E2eDL/FrameProducers/interface/predict_tf.h"
#include "tensorflow/core/graph/default_device.h"
#include "E2eDL/DataFormats/interface/FrameCollections.h"

e2e::Frame2D predict_tf(e2e::Frame4D& vinputFrame, string model_filename, string input_layer_name, string output_layer_name){
 tensorflow::Session* session;
 tensorflow::GraphDef graph_def;
 tensorflow::SessionOptions opts;
 std::vector<tensorflow::Tensor> outputs; // Store outputs
 // create a new session
 TF_CHECK_OK(NewSession(opts, &session));
 
 std::string graph_definition="E2eDL/"+model_filename;
 std::cout<<" >> Running Inference."<<endl;
 int batch_sz = vinputFrame.size();
 int frame_height = vinputFrame[0].size();
 int frame_width = vinputFrame[0][0].size();
 int no_channels = vinputFrame[0][0][0].size();
 //TF_CHECK_OK(ReadBinaryProto(Env::Default(), graph_definition, &graph_def));
 // load the graph definition, i.e. an object that contains the computational graph
 tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(graph_definition);
 tensorflow::Tensor tmp(tensorflow::DT_FLOAT, tensorflow::TensorShape({frame_height, frame_width}));
 /*auto _XTensor = tmp.matrix<float>();
 //std::copy_n(vEB_frame.begin(), vEB_frame.size(), tmp.flat<float>().data());
 for (int frame_row=0;frame_row<int(vinputFrame.size());frame_row++){
  for (int frame_col=0;frame_col<int(vinputFrame[0].size());frame_col++){
   _XTensor(frame_row,frame_col)=vinputFrame[frame_row][frame_col];
  }
 }*/
  
 tensorflow::Tensor x(tensorflow::DT_FLOAT, tensorflow::TensorShape({batch_sz, frame_height, frame_width, no_channels}));
 auto _XTensor = x.tensor<float,4>();
 for (int batch_idx=0; batch_idx<batch_sz; batch_idx++){
  for (int row_idx=0;row_idx<frame_height; row_idx++){
   for (int col_idx=0; col_idx<frame_width; col_idx++){
    for (int depth_idx=0; depth_idx<no_channels; depth_idx++){
     _XTensor(batch_idx, row_idx, col_idx, depth_idx) = vinputFrame[batch_idx][row_idx][col_idx][depth_idx];
    }
   }
  }
 }
  /*if(!x.CopyFrom(tmp, tensorflow::TensorShape({1, frame_height, frame_width, 1}))){
    std::cout<<" >> Reshape not successfull."<<endl;
  }*/
 // Set GPU options
 //graph::SetDefaultDevice("/gpu:0", &graph_def);
 //opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);
 //opts.config.mutable_gpu_options()->set_allow_growth(true);
 
 //int GPUID = std::stoi(params->getGpuDeviceStr());
 //setenv("CUDA_VISIBLE_DEVICES", "", GPUID);

 //std::cout << "Initial  visible_device_list : "<<opts.config.gpu_options().visible_device_list() << std::endl;
 opts.config.mutable_gpu_options()->set_allow_growth(true);
 opts.config.mutable_gpu_options()->set_per_process_gpu_memory_fraction(0.5);//params->getGpuMemoryRatio());
 
 
 // Load graph into session
 //TF_CHECK_OK(session->Create(graph_def));
 
 // create a session
 session = tensorflow::createSession(graphDef);
 
 // Initialize our variables
 
 TF_CHECK_OK(session->Run({{input_layer_name/*"inputs"*/, x}/*, {"y", y}*/}, {output_layer_name/*"softmax_1/Sigmoid"*/}, {}, &outputs)); // Get output
 //tensorflow::run(session, { { "x", x }, {"y", y} }, { "cost" }, &outputs);
 
 
 std::cout<<" >> Batch size is: "<<outputs.size()<<std::endl;
 e2e::Frame2D tmp_outputs = {{0,0},{0,0}};
 outputs.clear();
  
 session->Close();
 delete session;
 std::cout<<" >> Classification done"<<endl;
 // cleanup
 //tensorflow::closeSession(session);
 //delete graphDef;
 /*if (classifier_out>0.5){return 1;}
 else {return 0;}*/
 return tmp_outputs;
}
