//+
SetFactory("OpenCASCADE");
Merge "halfmodel.stp";
Recursive Delete {
  Surface{41}; 
}
//+
Delete {
  Surface{41}; 
}
//+
Delete {
  Surface{14}; Surface{41}; 
}

//+
Recursive Delete {
  Surface{7}; 
}
//+
Split Curve(40) {};
