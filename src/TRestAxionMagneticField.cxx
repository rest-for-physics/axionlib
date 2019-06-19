/*

// Don't forget to add the path of RestAxionLib as an environment varialbe (RestAxionLib_PATH)

<TRestAxionMagneticField> // cest comme le TRestReadout
// je pourrais tres bien def les parametres ici evidemment
	<addMagneticVolume type= "string" fileName="string" position="(posX,posY,posZ)" rotation="(rotX,rotY,rotZ)" />
<TRestAxionMagneticField/>

For now we consider that the bounds of the Volume are X,Y: -0.24 to 0.24 and Z:-5 to 5 
	
- type : 
	This is the type of mesh you want to choose for the geometry and for the regularity or not of the mesh.
	For now, there is the only "RegParallelepiped" option which is for a parallelepiped geometry of the magnetic Volume 
and a regular mesh. Maybe later there will be irregular mesh or other geometrical volume.

- fileName : this is the file from where the values of the magnetic field are read. The files are in data/magneticField. They all contain 3 columns for the position in the Volume and 3 columns to define the magnetic field vector. It is not the full path but only the name of the file.

- position, rotation : 
	By convention, the volume is build at the center of the frame. If you want to move it by translation, choose none zero posX, posY and/or posZ. If you want to move it by rotation, do the same for rotX, rotY, rotZ. rotX, rotY, rotZ correspond to Euler angles according to X,Y,Z (with the cannonical writing rotX->\psi, rotY->\phi, rotZ->\theta).
	//!\\ Be careful : for now, let the position and rotation on 0 because it is not still implemented (problem : no rotation or translation functions to my knowledge --> Create a temporal file which modifies the values of X Y Z (three first columns) of the fileName ??) 

// Later, maybe add a section to visualize the magnetic volume.	

*/

#include "TRestAxionMagneticField.h"
#include "TVectorD.h"
#include <iostream>

using namespace std;


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
#endif
}

void TRestAxionMagneticField::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

void TRestAxionMagneticField::LoadMagneticVolumeRegPar() {
#if defined USE_Garfield
	Int_t i=14;
	string sizemesh_s;
	sizemesh_s=fFileName[i];
	sizemesh_s=sizemesh_s+".";
	while (i<=18) {
		sizemesh_s=sizemesh_s+fFileName[i];
		i=i+1;}
	Double_t sizemesh=stod(sizemesh_s);
	TVectorD pos(3);  
	pos[0]=fPosX;
	pos[1]=fPosY;
	pos[2]=fPosZ;
	TVectorD rot(3);    
	rot[0]=fRotX;
	rot[1]=fRotY;
	rot[2]=fRotZ;
	const unsigned int nx=2*fXmax/sizemesh+1;
	const unsigned int ny=2*fYmax/sizemesh+1;
	const unsigned int nz=2*fZmax/sizemesh+1;
	Garfield::ComponentVoxel *fMesh=new ComponentVoxel();	
	//mesh->SetMesh(nx,ny,nz,fXmin,fXmax,fYmin,fYmax,fZmin,fZmax); // without any rotation and position consideration
	fMesh->SetMesh(nx,ny,nz,(fXmin+pos[0])*(cos(rot[1])*cos(rot[0])-sin(rot[1])*cos(rot[2])*sin(rot[0])+cos(rot[1])*sin(rot[0])+sin(rot[1])*cos(rot[2])*cos(rot[0])+sin(rot[1])*sin(rot[2])),(fXmax+pos[0])*(cos(rot[1])*cos(rot[0])-sin(rot[1])*cos(rot[2])*sin(rot[0])+cos(rot[1])*sin(rot[0])+sin(rot[1])*cos(rot[2])*cos(rot[0])+sin(rot[1])*sin(rot[2])),(fYmin+pos[1])*(-sin(rot[1])*cos(rot[0])-cos(rot[1])*cos(rot[2])*sin(rot[0])-sin(rot[1])*sin(rot[0])+cos(rot[1])*cos(rot[2])*cos(rot[0])+cos(rot[1])*sin(rot[2])),(fYmax+pos[1])*(-sin(rot[1])*cos(rot[0])-cos(rot[1])*cos(rot[2])*sin(rot[0])-sin(rot[1])*sin(rot[0])+cos(rot[1])*cos(rot[2])*cos(rot[0])+cos(rot[1])*sin(rot[2])),(fZmin+pos[2])*(sin(rot[2])*sin(rot[0])-sin(rot[2])*cos(rot[0])+cos(rot[2])),(fZmax+pos[2])*(sin(rot[2])*sin(rot[0])-sin(rot[2])*cos(rot[0])+cos(rot[2]))); // to future translation and rotation. Let rota and pos on 0 for now.	
	fMesh->LoadMagneticField((TString)getenv("RestAxionLib_PATH")+"data/magneticField"+fFileName,"XYZ",1,1);
	fMesh->EnableInterpolation(true);
#else
    cout << "This REST is not complied with garfield, it cannot load any magnetic field Volume!" << endl;
#endif
}



void TRestAxionMagneticField::InitFromConfigFile() {
    this->Initialize();
    string gasbVolume;
    size_t position = 0;
    while ((gasbVolume = GetKEYStructure("addMagneticVolume",position)) != "") { // debugs or warnings later
	fType=GetParameter("type","RegParallelepiped");
	if(fType=="RegParallelepiped")
		fFileName=GetParameter("fileName","magneticField_004_Bykovskiy_201906.dat");
		fPosX=StringToDouble(GetParameter("posX","0"));
		fPosY=StringToDouble(GetParameter("posY","0"));
		fPosZ=StringToDouble(GetParameter("posZ","0"));
		fRotX=StringToDouble(GetParameter("rotX","0"));
		fRotY=StringToDouble(GetParameter("rotY","0"));
		fRotZ=StringToDouble(GetParameter("rotZ","0"));
		LoadMagneticVolumeRegPar(); 
	//else{} // To see later
	}
}


void TRestAxionMagneticField::PrintMetadata() {
     TRestMetadata::PrintMetadata();
     metadata << " - Type of mesh : " << fType << endl;
     metadata << " - File loaded : " << fFileName<< endl;
     metadata << " - Position : " << "(" << fPosX << "," << fPosY << "," << fPosZ << ")"<< endl;
     metadata << " - Rotation : " << "(" << fRotX << "," << fRotY << "," << fRotZ << ")"<< endl;
     metadata << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
