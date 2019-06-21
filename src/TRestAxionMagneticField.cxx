/*

// Don't forget to add the path of RestAxionLib as an environment varialbe (RestAxionLib_PATH)

<TRestAxionMagneticField> 
	<addMagneticVolume fileName="string" position="fPos" />
<TRestAxionMagneticField/>

For now we consider that the bounds of the Volume are X,Y: -0.24 to 0.24 and Z:-5 to 5 

- fileName : this is the file from where the values of the magnetic field are read. The files are in data/magneticField. They all contain 3 columns for the position in the Volume and 3 columns to define the magnetic field vector. It is not the full path but only the name of the file.

- position : 
	By convention, the volume is build at the center of the frame. If you want to move it by translation, choose put a non null vector position "fPos".

// Later, maybe add a section to visualize the magnetic volume (Can be done with ViewField())

// We try to do rotations but for now we could not so we choose to define the reference from the center of the magnetic bore configuration and the rotations will be done by the solar classes. 

*/

#include "TRestAxionMagneticField.h"

using namespace std;
using namespace Garfield;


ClassImp(TRestAxionMagneticField);


TRestAxionMagneticField::TRestAxionMagneticField() : TRestMetadata() { Initialize(); }


TRestAxionMagneticField::TRestAxionMagneticField(const char* cfgFileName, string name) : TRestMetadata(cfgFileName) {
    cout << "Entering TRestAxionSolarModel constructor( cfgFileName, name )" << endl;    
     
    Initialize();

    LoadConfigFromFile(fConfigFileName, name);

    PrintMetadata();
}


TRestAxionMagneticField::~TRestAxionMagneticField() {
    //TRestAxionMagneticField destructor

    debug << "Entering ... TRestAxionMagneticField() destructor." << endl;

#if defined USE_Garfield
    delete fMesh;
    delete fSetOfField ;
#endif
}

void TRestAxionMagneticField::Initialize() { 
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
    
    fFileName="magneticField_Bykovskiy_201906.dat";
    fNofVolumes=0;
    sizeMesh=0.04;
    fXmin=-0.24;
    fXmax=0.24; 
    fYmin=-0.24;
    fYmax=0.24; 
    fZmin=-5.0;
    fZmax=5.0; 
    for(int i=0;i<3;i++) fPos[i]=0.0;
    // How to initialize fPos and fPositions ?

#if defined USE_Garfield
    fSetOfField = new Garfield::Sensor();
    fMesh = new Garfield::ComponentVoxel();
#else
    fSetOfField = NULL;
    fMesh = NULL;
#endif
}

void TRestAxionMagneticField::LoadMagneticVolume(TVector3 Pos) {

#ifdef USE_Garfield

	/*** Read information from the fileName ***/

	fstream file_r; 
	file_r.open((string)getenv("RestAxionLib_PATH")+"/data/magneticField/"+(string)fFileName); // original file 
	if (file_r.is_open()) 
	{
	  ofstream file_w;
    	  file_w.open("/tmp/tmp_bField.txt"); // temporal file for the translation and without the first line
	  TVector3 position; // translation vector
	  double read;
	  int i=0;
	  int j=0;
	  while (file_r >> read) 
    	  {  
		  if(i<4) {
			  if (i==0) { double fXmax=read; double fXmin=-fXmax; i=i+1;}
		 	  else if (i==1) { double fYmax=read; double fYmin=-fYmax; i=i+1;}
			  else if (i==2) { double fZmax=read; double fZmin=-fZmax; i=i+1;}
			  else if (i==3) { double sizeMesh=read; i=i+1;}
			
			  }		
		
		  else {
			  if(j<3) {
				  position[j]=read;
				  if(j==2) {
				  	  for(int k=0;k<3;k++) {
						  position[k]=position[k]+Pos[k];
						  file_w<<position[k];
						  file_w<<"\t";}
					   }
				  j=j+1;
				  }
		
			  else {
				  if(j==5) {
					  file_w<<read;
					  file_w<<"\n";
					  j=0;
					   }
				  else {
					  file_w<<read;
					  file_w<<"\t";
					  j=j+1;
				       }
			        }
		       }				
				     
    	  } 
	file_w.close();
    	file_r.close();

	/*** Create the mesh if file is open ***/
	int nx=(int)(2*fXmax/sizeMesh)+1;
	int ny=(int)(2*fYmax/sizeMesh)+1;
	int nz=(int)(2*fZmax/sizeMesh)+1;	
	//cout << fXmin+Pos[0]<<" " <<fXmax+Pos[0]<< " " << fYmin+Pos[1] << " " << fYmax+Pos[1] << " " << fZmin+Pos[2]<< " " <<fZmax+Pos[2]<<endl;
	fMesh->SetMesh(nx,ny,nz,fXmin+Pos[0],fXmax+Pos[0],fYmin+Pos[1],fYmax+Pos[1],fZmin+Pos[2],fZmax+Pos[2]);

	/*** Fill the mesh with the magnetic field if file is open ***/
	fMesh->LoadMagneticField("/tmp/tmp_bField.txt","XYZ",1,1); // pour le nom a changer plus tard pour eviter  d effacer a chaque fois
	fMesh->EnableInterpolation(true);
	}
	else cout <<" Cannot find the file "<< endl;
    	

	
#else
    cout << "This REST is not complied with garfield, it cannot load any magnetic field Volume!" << endl;
#endif

}


void TRestAxionMagneticField::AddFieldComponent() { 
	fSetOfField->AddComponent(fMesh);
	fNofVolumes ++;
}


void TRestAxionMagneticField::InitFromConfigFile() {
    this->Initialize();
    string bVolume;
    size_t position = 0;
    while ((bVolume = GetKEYDefinition("addMagneticVolume",position)) != "") // debugs or warnings later
    { 
	fFileName=GetParameter("fileName","magneticField_Bykovskiy_201906.dat");
	fPos=StringTo3DVector(GetParameter("position","(0,0,0)"));
	LoadMagneticVolume(fPos);
	AddFieldComponent(); 
	fPositions.push_back(fPos);
    }
}


void TRestAxionMagneticField::PrintMetadata() {
     
     TRestMetadata::PrintMetadata();

     metadata << " - File loaded : " << fFileName<< endl;
     metadata << " - Size of the mesh :" << sizeMesh << endl;
     metadata << " - Bounds : " << endl;
     metadata << "    * xmin = " << fXmin << ", xmax = " << fXmax << endl;        
     metadata << "    * ymin = " << fYmin << ", ymax = " << fYmax << endl;
     metadata << "    * zmin = " << fZmin << ", zmax = " << fZmax << endl;
     metadata << " - Number of magnetic volumes : "<< fNofVolumes << endl; 
     metadata << " ------------------------------------------------ " << endl;
     metadata << " - Positions of Volumes :" << endl;
     int n=0;
     double x,y,z;
     for (int p = 0; p <fNofVolumes; p++) {
	x=fPositions[p][0];
	y=fPositions[p][1];
	z=fPositions[p][2];
	metadata << "    * Volume " << p+1 << " : " << "("<<x<<","<<y<<","<<z<<")" << endl;} // Bcz I don't know if the order will count later
     metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

}
