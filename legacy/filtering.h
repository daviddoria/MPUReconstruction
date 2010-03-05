
extern std::string IntToString(int i);
extern void extractSlice( vtkPoints *, vtkIntArray *, vtkPoints *, int);
extern void VolumeInit( vtkImageData *, vtkPoints *, vtkIntArray *);
extern void ImageInit(vtkImageData *, int *, double *);
extern void BinarizeVolume(vtkImageData *, vtkPoints *, vtkIntArray *);
extern void BinarizeSlice(vtkImageData *, vtkPoints *, vtkFloatArray *,int);
extern void GaussianSmooth(vtkImageData *);
extern void NormalEstimation(vtkImageData *);
extern void ExtractNormal(vtkImageData *, vtkPoints *, vtkFloatArray *);
extern void AssociatePointNormal(vtkPoints *, vtkFloatArray *, vtkPolyData*);
