#ifndef IO_H
#define IO_H
// extern void Load(vtkPoints *, vtkIntArray *, const char*);
// extern void ImageWrite( vtkImageData *, const char *);
// extern void PolyWrite( vtkPolyData *, const char *);
// extern void PWNWrite( vtkPoints *, vtkFloatArray *, const char*);
// 
void Load(vtkPoints *, vtkIntArray *, const char*);
void ImageWrite( vtkImageData *, const char *);
void PolyWrite( vtkPolyData *, const char *);
void PWNWrite( vtkPoints *, vtkFloatArray *, const char*);
void MPULoad( PointSet *ps, const char *file);

#endif