# AOGDetector
                          Implementation of the papers

        "Discriminatively Trained And-Or Tree Model for Object Detection"
        CVPR, 2013.        

        "Modeling Occlusion by Discriminative AND-OR Structures"
        ICCV, 2013        

        "Integrating Context and Occlusion for Car Detection by Hierarchical And-Or Model"
        ECCV, 2014

Copyright (c) 2014 Center for Vision, Cognition, Learning and Autonomy at UCLA <http://vcla.stat.ucla.edu/>
            Written by Matt Tianfu Wu (primary contact) <tfwu@stat.ucla.edu> 
                       Xi Song                          <xsong.lhi@gmail.com> (object detection module)                       
                       Bo Li                            <boli.lhi@gmail.com>  (occlusion module)                     


                       Code Release Log

10-14-2014: Release detection code  (full code will be released in future)                    


                       MODULE 1: object detection 

      VOC Configuration
  All the configuration can be specified in voc_config.xml in "<PATH_TO_RGM_Source_Code>/config" folder.

      VOC Testing
  Use Terminal and cd to the <PATH_TO_Build>/bin where the executable "Entry" is built
  Run:   Entry vocEval <PATH_To_voc_config_xml>
  e.g.,  Entry vocEval /home/tfwu/Code/RGM/config/voc_config.xml

     DetectionInASingleImage
  Use Terminal and cd to the <PATH_TO_Build>/bin where the executable "Entry" is built
  Run:   Entry detect <PATH_To_voc_config_xml> 
  e.g.,  Entry detect /home/tfwu/Code/RGM/config/voc_config.xml

     BatchDetectionInDirectory
  Use Terminal and cd to the <PATH_TO_Build>/bin where the executable "Entry" is built
  Run:   Entry detectBatch <PATH_To_voc_config_xml> 
  e.g.,  Entry detectBatch /home/tfwu/Code/RGM/config/voc_config.xml  

     ConvertVOCRel5Model
  Use Terminal and cd to the <PATH_TO_Build>/bin where the executable "Entry" is built
  Run:   Entry cvtVOCRel5Model <PATH_To_VOCRel5_Model_Mat_File> 
  e.g.,  Entry cvtVOCRel5Model /home/tfwu/Code/RGM/externalModels/VOC2007/bicycle_final.mat
  The converted model will be saved (e.g., /home/tfwu/Code/RGM/data/bicycle_final.mat.bin),
  which you can used to test detection.
  Note: before using this to convert a matlab model trained using VOC-Release 5 from Dr. Ross Girshick, 
        please use <PATH_TO_RGM_Source_Code>/matlab/convertVOCR5.m to add a "DFS" field to the model structure.
