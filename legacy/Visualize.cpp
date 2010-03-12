
#include "main.h"
#include "Visualize.h"
#include "vtkDataSetMapper.h"
#include "vtkDataSet.h"

void PointsActor( vtkPoints *points, vtkActor ** ppointsActor, double color[3])
{
	// Create vtkCellArray
	vtkCellArray *verts = vtkCellArray::New();
	for( int i=0; i<points->GetNumberOfPoints(); i++)
		verts->InsertNextCell(1, &i );

	// Create a PolyData dummy variable
	vtkPolyData *dummyPolyData = vtkPolyData::New(); 
	dummyPolyData->SetPoints( points );
	dummyPolyData->SetVerts( verts );

	vtkPolyDataMapper *pointsMapper = vtkPolyDataMapper::New();
	pointsMapper->SetInput( dummyPolyData );
	pointsMapper->ScalarVisibilityOff();
	(*ppointsActor)->SetMapper( pointsMapper );
	(*ppointsActor)->GetProperty()->SetColor( color );

	pointsMapper->Delete();
	dummyPolyData->Delete();
	verts->Delete();
}

void RawImageActor( vtkImageData *image, vtkActor2D **pimageActor)
{
	float *pValue = (float *)image->GetScalarPointer(); 
	int offset;
	float maxValue = 0;
	float minValue = 1000;	
	for (int i=0; i<512;i++)
		for(int j=0; j<512; j++){
			offset = j*512+i;
			if (pValue[offset]>maxValue)
				maxValue = pValue[offset];
			if (pValue[offset]<minValue)
				minValue = pValue[offset];
		}
	int colorWindow = maxValue-minValue;
	int colorLevel = (maxValue+minValue)/2;
	vtkImageMapper * imageMapper = vtkImageMapper::New();
	imageMapper->SetInput( image );
	imageMapper->SetColorLevel( colorLevel );
	imageMapper->SetColorWindow( colorWindow );

	(*pimageActor)->SetMapper( imageMapper );	

}

void PolyActor( vtkPolyData * poly, vtkActor **ppolyActor, double color[3])
{
	vtkPolyDataMapper *pointsMapper = vtkPolyDataMapper::New();
	pointsMapper->SetInput( poly );
	pointsMapper->ScalarVisibilityOff();
	(*ppolyActor)->SetMapper( pointsMapper );
	(*ppolyActor)->GetProperty()->SetColor( color );
	pointsMapper->Delete();
}

void Render(vtkRenderer **pren)
{
	// Renderer and RenderWindow
	(*pren)->SetBackground( 0, 0, 0);
	(*pren)->ResetCamera();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->SetSize( 512,512 );
	renWin->AddRenderer( *pren );
	vtkRenderWindowInteractor *iren =	vtkRenderWindowInteractor::New();
	iren->SetRenderWindow( renWin );

	renWin->Render();

	iren->Initialize();
	iren->Start();

	iren->Delete();
	renWin->Delete();
}

void PointsVisualize(vtkPoints *points1)
{
	vtkActor *points1Actor = vtkActor::New();
	double color1[3] = {1,1,1};
	PointsActor( points1, &points1Actor, color1);
	vtkRenderer *ren = vtkRenderer::New();
	ren->AddActor(points1Actor);
	Render( &ren );
	ren->Delete();
	points1Actor->Delete();
}

void Points2Visualize(vtkPoints *points1, vtkPoints *points2)
{
	vtkActor *points1Actor = vtkActor::New();
	vtkActor *points2Actor = vtkActor::New();
	double color1[3] = {1,0,0};
	double color2[3] = {0,1,0};
	PointsActor( points1, &points1Actor, color1);
	PointsActor( points2, &points2Actor, color2);
	vtkRenderer *ren = vtkRenderer::New();
	ren->AddActor(points1Actor);
	ren->AddActor(points2Actor);
	Render( &ren );
	ren->Delete();
	points1Actor->Delete();
	points2Actor->Delete();
}

void ImageVisualize(vtkImageData *img)
{
//	vtkImageActor *imgActor = vtkImageActor::New();
//	imgActor->SetInput(m_WLFilter->GetOutput());
	vtkActor2D *imgActor = vtkActor2D::New();
	RawImageActor(img, &imgActor);
	vtkRenderer *ren = vtkRenderer::New();
	ren->AddActor( imgActor );
	Render( &ren );
	ren->Delete();
	imgActor->Delete();
}

void VolumeVisualize(vtkImageData *img)
{
	vtkVolumeRayCastCompositeFunction *rayCast = vtkVolumeRayCastCompositeFunction::New();

	vtkVolumeRayCastMapper *volumeMapper = vtkVolumeRayCastMapper::New();
	volumeMapper->SetInput(img);
	volumeMapper->SetVolumeRayCastFunction(rayCast);

	vtkVolume *volume = vtkVolume::New();
	volume->SetMapper(volumeMapper);

	vtkRenderer *ren = vtkRenderer::New();
	ren->AddViewProp(volume);

	Render(&ren);

	rayCast->Delete();
	volumeMapper->Delete();
	volume->Delete();
	ren->Delete();
}