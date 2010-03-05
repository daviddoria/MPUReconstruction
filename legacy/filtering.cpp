#include "main.h"
#include "Visualize.h"
#include "IO.h"

std::string IntToString(int i)
{
	std::ostringstream sstr;
	sstr<<i;	
	std::string s = sstr.str();
	return s;
}

void extractSlice(vtkPoints *points, vtkIntArray *idxPts, vtkPoints *pointsSlice, int sliceIdx)
{
	int start=0, end=0;
	if (sliceIdx >0)
		for (int i=0; i<sliceIdx; i++){
			start =  start + idxPts->GetComponent(1, i);
		}
	end = start + idxPts->GetComponent(1, sliceIdx);
	vtkIdList *ls = vtkIdList::New();
	for (int i=start; i<end; i++)
		ls->InsertNextId(i);

	points->GetPoints(ls, pointsSlice);
	ls->Delete();
}


void VolumeInit( vtkImageData *contourImg, vtkPoints *points, vtkIntArray *idxPts)
{
	double *bounds = points->GetBounds();
	int xsize = vtkMath::Round(bounds[1])-vtkMath::Round(bounds[0])+40;
	int ysize = vtkMath::Round(bounds[3])-vtkMath::Round(bounds[2])+40;
	int zsize = idxPts->GetNumberOfComponents()+10;
	int xspacing = 1;
	int yspacing = 1;
	int zspacing = (vtkMath::Round(bounds[5])-vtkMath::Round(bounds[4]))/(idxPts->GetNumberOfComponents()-1);
	contourImg->SetDimensions(xsize,ysize,zsize);
	contourImg->SetOrigin(vtkMath::Round(bounds[0])-20,vtkMath::Round(bounds[2])-20,vtkMath::Round(bounds[4])-5*zspacing);
	contourImg->SetSpacing(xspacing,yspacing,zspacing);
	contourImg->SetScalarTypeToFloat();
	contourImg->AllocateScalars();
}

void ImageInit(vtkImageData *img, int *dim, double *origin)
{
	(img)->SetOrigin(origin);
	(img)->SetDimensions(dim);
	(img)->SetSpacing(1,1,1);
	(img)->SetScalarTypeToFloat();
	(img)->AllocateScalars();
}

void BinarizeSlice( vtkImageData *contourImg,vtkPoints *pointsSlice,vtkFloatArray *scalars,int z)
{
	int *dim = contourImg->GetDimensions();
	double *origins = contourImg->GetOrigin();
	int dimSlice[3]; dimSlice[0]=dim[0]; dimSlice[1]=dim[1]; dimSlice[2]=1;
	double originSlice[3]; originSlice[0]=origins[0]; originSlice[1]=origins[1]; originSlice[2]=0;
	vtkImageData *imgSlice = vtkImageData::New();
	ImageInit(imgSlice,dimSlice,originSlice);
	vtkFloatArray *scalarSlice = vtkFloatArray::New(); scalarSlice->SetNumberOfTuples(dim[0]*dim[1]);
	// convert pointsSlice to double array 
	int numPts = pointsSlice->GetNumberOfPoints()-1; // do not repeat the last point
	double *pts_double = new double [numPts*3];
	vtkPoints *pt = vtkPoints::New();
	for (int i=0; i< numPts ; i++)
	{
		double *p = pointsSlice->GetPoint( i );
		pts_double[i*3+0] = p[0];
		pts_double[i*3+1] = p[1];
		pts_double[i*3+2] = p[2];
	}	
	double z_value = vtkMath::Round(pts_double[2]);
	// pointinpolygon calculation, set corresponding scalars region
	int offset;
	double n[3];
	double *b = pointsSlice->GetBounds();
	b[4]=z_value; b[5]=z_value;
	vtkPolygon::ComputeNormal(pointsSlice, n);
	int in_poly; // flag to indicate whether in/out polygon
	vtkPoints *pin = vtkPoints::New();
	for (int j=0; j<dim[1]; j++)
		for (int i=0; i<dim[0]; i++)
		{
			double x[3] = {i+origins[0], j+origins[1], z_value};
			offset = i+j*dim[0]+(z+5)*dim[0]*dim[1];			
			in_poly = vtkPolygon::PointInPolygon( x, numPts, pts_double, b, n);
//			cout<<offset<<"   "<<in_poly<<endl;
			if (in_poly){
				scalars->SetComponent(offset,0,100);				
				pin->InsertNextPoint(x[0],x[1],x[2]);
				scalarSlice->SetComponent(i+j*dim[0],0,100);
//				Points2Visualize(pointsSlice,pt);
			}
			else{
				scalars->SetComponent(offset,0,0);
				scalarSlice->SetComponent(i+j*dim[0],0,0);
			}
//			pt->InsertNextPoint(x[0],x[1],x[2]);			
		}
//	Points2Visualize( pointsSlice,pin);
	imgSlice->GetPointData()->SetScalars(scalarSlice);		
	std::string	filename = "Binary" + IntToString(z) + ".vti";
	//ImageWrite(imgSlice, filename.c_str());
	pin->Delete();
	pt->Delete();
	delete []pts_double;
	scalarSlice->Delete();
	imgSlice->Delete();
}

void BinarizeVolume(vtkImageData *contourImg, vtkPoints * points, vtkIntArray *idxPts)
{
	VolumeInit(contourImg, points, idxPts); // dim, origin, spacing, ...
	int offset;	
	int numSlice = idxPts->GetNumberOfComponents();
	int *dim = contourImg->GetDimensions();
	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->SetNumberOfTuples(dim[0]*dim[1]*dim[2]);
	for(int i=0; i<dim[0]; i++)
		for(int j=0; j<dim[1]; j++)
			for(int k=0; k<dim[2]; k++){
				offset = i+ j*dim[0] + k*dim[0]*dim[1];
				scalars->SetComponent(offset,0,0);
			}
	
	for (int i=0; i<numSlice; i++){
		vtkPoints *pointsSlice = vtkPoints::New();  // ??????????????????????
		extractSlice(points,idxPts,pointsSlice,i); 
//		PointsVisualize( pointsSlice );
		BinarizeSlice(contourImg, pointsSlice,scalars, i); 
		pointsSlice->Delete();
	}	
	contourImg->GetPointData()->SetScalars(scalars);
	scalars->Delete();
}

void GaussianSmooth(vtkImageData *contourImg)
{
	vtkImageGaussianSmooth *gaussianFilter = vtkImageGaussianSmooth::New();
	gaussianFilter->SetInput(contourImg);
	gaussianFilter->SetDimensionality(3);
	gaussianFilter->SetStandardDeviations(3,3,15);
//	gaussianFilter->SetRadiusFactor(15);
	gaussianFilter->Update();
	(contourImg)->Initialize();
	(contourImg)->DeepCopy(gaussianFilter->GetOutput());
	gaussianFilter->Delete();
}

void NormalEstimation(vtkImageData *smoothImg)
{
	vtkImageGradient *gradientFilter = vtkImageGradient::New();
	gradientFilter->SetInput(smoothImg);
	gradientFilter->SetDimensionality(3);
	gradientFilter->HandleBoundariesOn();
	gradientFilter->Update();
	(smoothImg)->Initialize();
	(smoothImg)->DeepCopy(gradientFilter->GetOutput());
	gradientFilter->Delete();
}

void ExtractNormal(vtkImageData *gradImg, vtkPoints *points, vtkFloatArray *normal)
{
	int numPts = points->GetNumberOfPoints();
	int *dim = gradImg->GetDimensions();
	double *spacing = gradImg->GetSpacing();
	double *origin = gradImg->GetOrigin();
	(normal)->SetNumberOfComponents(3);
	(normal)->SetNumberOfTuples(numPts);
	vtkDataArray *grad = gradImg->GetPointData()->GetScalars();	
	int offset;
	double g[3];
	double *p;
	double mag;
	for (int i=0; i<numPts; i++){
		p=points->GetPoint(i);
		p[0]=(vtkMath::Round(p[0]-origin[0]))/spacing[0]; 
		p[1]=(vtkMath::Round(p[1]-origin[1]))/spacing[1]; 
		p[2]=(vtkMath::Round(p[2]-origin[2]))/spacing[2];
		offset = p[0]+p[1]*dim[0]+p[2]*dim[0]*dim[1];
		g[0]=-grad->GetComponent(offset,0); 
		g[1]=-grad->GetComponent(offset,1); 
		g[2]=-grad->GetComponent(offset,2);
		mag = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
		g[0]=g[0]/mag;		g[1]=g[1]/mag;      g[2]=g[2]/mag;
		(normal)->SetComponent(i,0,g[0]);
		(normal)->SetComponent(i,1,g[1]);
		(normal)->SetComponent(i,2,g[2]);
	}
}

void AssociatePointNormal(vtkPoints * points, vtkFloatArray * normal, vtkPolyData* poly)
{
	(poly)->SetPoints( points );
	(poly)->GetPointData()->SetNormals(normal);
	vtkCellArray *verts = vtkCellArray::New();
	for( int i=0; i<points->GetNumberOfPoints(); i++)
		verts->InsertNextCell(1, &i );
	(poly)->SetVerts(verts);
}
