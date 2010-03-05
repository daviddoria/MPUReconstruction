#include "main.h"
#include "filtering.h"
#include "IO.h"

void Load( vtkPoints *points, vtkIntArray *numPts, const char *fileName)
{
	ifstream DataFile;
	DataFile.open(fileName);
	if (!DataFile.is_open())
  {
		std::cerr<< "open error!"<<endl;
  }

	DataFile.seekg(0);
	char strTemp[28];
	char charTemp;
	double xcor, ycor, zcor;
	double x,y,z;
	int numSlice, numPoints;
	DataFile.get(strTemp,19); // Number of Slices
	DataFile.get(strTemp,3);    // Slices;
	numSlice = atoi(strTemp);
	DataFile.getline(strTemp,28);

	numPts->SetNumberOfComponents(numSlice);
	vtkPoints *poriginal = vtkPoints::New();
	for (int i=0; i<numSlice; i++){
		DataFile.getline(strTemp,28);
		DataFile.get(strTemp,12);
		DataFile.get(strTemp,5);
		numPoints = atoi(strTemp);
		for (int j=0; j<numPoints;j++){
			DataFile.get(strTemp,28,',');
			x = atof(strTemp);
			DataFile.get(charTemp);
			DataFile.get(strTemp,28,',');
			y = atof(strTemp);			
			DataFile.get(charTemp);
			DataFile.get(strTemp,28,'\n');
			z = atof(strTemp);
			xcor = x*10+256; ycor = y*10+256; zcor = z*2; //compensation
			points->InsertNextPoint(xcor,ycor,zcor);
			(poriginal)->InsertNextPoint(x,y,z);
		}
		numPts->InsertComponent(1,i,numPoints);
		DataFile.getline(strTemp,28);
	}
	DataFile.close();
	vtkPolyData *poly = vtkPolyData::New();

}

void MPULoad( PointSet *ps, const char *file)
{
	ifstream DataFile;
	DataFile.open(file);
	if (!DataFile.is_open())
		std::cerr<< "open error!"<<endl;

	DataFile.seekg(0);
	char strTemp[28];
	char charTemp;
	double x,y,z;	
	int numPoints;
	DataFile.get(strTemp,8);
	numPoints = atoi(strTemp);
	ps->setPointSize(numPoints); //numPts

	////read points
	for(int i=0; i<numPoints; i++){
		DataFile.get(strTemp,28,' ');
		x = atof(strTemp);
		DataFile.get(charTemp);
		DataFile.get(strTemp,28,' ');
		y = atof(strTemp);			
		DataFile.get(charTemp);
		DataFile.get(strTemp,28,'\n');
		z = atof(strTemp);
		ps->setPoint(i, x, y, z);
	}
	////read normals
	for(int i=0; i<numPoints; i++){
		DataFile.get(strTemp,28,' ');
		x = atof(strTemp);
		DataFile.get(charTemp);
		DataFile.get(strTemp,28,' ');
		y = atof(strTemp);			
		DataFile.get(charTemp);
		DataFile.get(strTemp,28,'\n');
		z = atof(strTemp);
		ps->setNormal(i, x, y, z);
	}

	DataFile.close();	
}

std::string DoubleToString(double i)
{
	std::ostringstream sstr;
	sstr<<i;	
	std::string s = sstr.str();
	return s;
}

void PWNWrite( vtkPoints *points, vtkFloatArray *normal, const char* fileName)
{
	ofstream DataFile;
	DataFile.open(fileName);	
	int numPts = points->GetNumberOfPoints();
	std::string s = DoubleToString(numPts);
	DataFile<<s;		   //write # of points
	DataFile<<'\n';
	double p[3];
	double n[3];
	for (int i=0; i<numPts;i++) //write points first
	{
		points->GetPoint(i,p);
		s = DoubleToString((p[0]-256)/10);
		DataFile<<s;
		DataFile<<' ';
		s = DoubleToString((p[1]-256)/10);
		DataFile<<s;				
		DataFile<<' ';
		s = DoubleToString(p[2]/10);
		DataFile<<s;
		DataFile<<'\n';
	}
	for (int i=0; i<numPts;i++) //write normals second
	{
		n[0]=normal->GetComponent(i,0);
		n[1]=normal->GetComponent(i,1);
		n[2]=normal->GetComponent(i,2);
		s = DoubleToString(n[0]);
		DataFile<<s;
		DataFile<<' ';
		s = DoubleToString(n[1]);
		DataFile<<s;				
		DataFile<<' ';
		s = DoubleToString(n[2]);
		DataFile<<s;
		DataFile<<'\n';
	}
}
