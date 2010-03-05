/*
	Contour Structure Frame:

	Structures.dcm --
		StructureSetROISequence --  3006 0020
			SQItem [n] --
				ROIName -- ( Match the ROI with the organs requested, return index n)  3006 0026
		ROIContourSequence --   3006 0039
			SQItem [n] --
				ContourSequence -- 3006 0040
					SQItem [m] -- The mth contour, they may be on the same slice.
						ContourData -- the coordinates	3006 0050
						NumberofContourPoints -- # of points on this slice	3006 0046
						ContourImageSequence -- 3006 0016
							SQItem [1]--
								ReferencedSOPInstanceID --  Match contour with image 0008 1155					

*/



#include "main.h"

///////////////////////////////////   Load Sequence items such as contour and ID from DCM files uing GDCM

extern void ArrayForm2D(std::string, uint32_t, vtkPoints **);

void ContourLoad( vtkPoints** pPoints )
{
	gdcm::File *f = new gdcm::File();
	// organ requested
	std::vector <std::string> organ(3);
	organ[0]="CTVajm";
	organ[1]="rectum";
	organ[2]="bladder";

	// organ index according to Structure Set ROI Sequence
	std::vector <int> organ_idx(3);

	// Structure dcm file
	std::string fileName = "G:/UCSD/data/Patient1/P1_CBCT1/P1_CBCT1_Structures.dcm" ;
	uint16_t group = 0x3006; 

	//f = gdcm::File::New( );   
	f->SetFileName( fileName );

	// Make sure the Data Element will be loaded
	f->AddForceLoadElement ( group, 0x0020);
	f->AddForceLoadElement ( group, 0x0039);
	f->AddForceLoadElement ( group, 0x0050);
	f->SetMaxSizeLoadEntry(0xffff);
	f->Load();  


	// Find the organ index
	//StructureSetROISequence
	gdcm::SeqEntry *dicom_seq = f->GetSeqEntry(group,0x0020);
	int TotalNumSQ = dicom_seq->GetNumberOfSQItems();
	if (dicom_seq!=0)
	{				
		int CurrentSQNum=0;
		// loop all the StructureSetROISequence items to compare ROINames with organs
		while (CurrentSQNum<TotalNumSQ)
		{
			gdcm::SQItem *sqi=dicom_seq->GetSQItem(CurrentSQNum);
			std::string ROI_name = sqi->GetValEntry(group, 0x0026)->GetValue();

			//StructROI = "bladder "; remove tha last space
			if (ROI_name.find(' ')<ROI_name.size())
				ROI_name.erase(ROI_name.find(' '));
			for (int i=0;i<3;i++)
			{
				if (ROI_name.compare(organ[i])==0)
				{
					//std::cout<<ROI_name<<" "<<CurrentSQNum<<std::endl;
					organ_idx[i]=CurrentSQNum;
				}
			}			
			CurrentSQNum++;
		}
	}

	// Find the contour;
	// ROIContourSequence;
	//std::ofstream o;
	dicom_seq = f->GetSeqEntry(group,0x0039);	
	if (dicom_seq!=0)
	{
		// Get each organ's contour and contour ID
		for (int i=0;i<1;i++)
		{
			//std::string filename = "D:/research work/UCSD/data/Patient1/P1_CBCT1/" + organ[i] + ".txt";
			//o.open(filename.c_str());
			gdcm::SQItem *sqi = dicom_seq->GetSQItem(organ_idx[i]);		
			// ContourSequence
			gdcm::SeqEntry *ContourSeq = sqi->GetSeqEntry(group,0x0040);
			for (int j=3; j<4;j++)//j<ContourSeq->GetNumberOfSQItems(); j++)
			{	
				
				gdcm::SQItem *temp = ContourSeq->GetSQItem(j);
				std::string contour_slice = temp->GetValEntry(group, 0x0050)->GetValue();
				uint32_t len = temp->GetValEntry(group, 0x0050)->GetLength();;
				
				// Load into vtkFloatArray				
				ArrayForm2D(contour_slice,len, pPoints);				
			
				//std::string ContourID = ContourSeq->GetSQItem(j)->GetSeqEntry(group,0x0016)->
				//	GetFirstSQItem()->GetEntryString(0x0008,0x1155);					
//				o.write((char*)ContourSlice->GetBinArea(), ContourSlice->GetLength());	
								
			}
			//o.close();
		}
	}

	delete f;

}