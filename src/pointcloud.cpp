#include "../headers/pointcloud.h"
#include <cmath>
#include <liblas/liblas.hpp>
#include "../utils/simple_math.h"

using namespace std;


void DigitalElevModel::init(int _nx, int _ny, float _dx, float _dy, float _x0, float _y0){
	nx   = _nx;			
	ny   = _ny;			
	dx   = _dx;			
	dy   = _dy;			
	x0   = _x0;			
	y0   = _y0;
	
	xmin = x0 - dx/2;
	ymin = x0 - dx/2;
	xmax = (dx*nx) + xmin; 
	ymax = (dy*ny) + ymin;
	
	elev.resize(nx*ny);
	
}

void DigitalElevModel::initFromBounds(int _xmin, int _xmax, float _ymin, float _ymax, float _dx, float _dy){
	xmin = _xmin;			
	xmax = _xmax;			
	ymin = _ymin;			
	ymax = _ymax;			
	dx   = _dx;			
	dy   = _dy;
	nx   =  (xmax- xmin)/ dx ; 
	ny   =  (ymax- ymin)/ dy ; 
	x0   =  xmin + dx/2; 
	y0   =  ymin + dy/2;
	
	elev.resize(nx*ny);
	 			
}

vector<int> DigitalElevModel::get_index(float x, float y){
	vector<int>index(2);
	index[0] = floor((x- xmin)/ dx);
	index[1] = floor((y- ymin)/ dy);
	return index;
}


void DigitalElevModel::printToFile(string filename){
	ofstream fout(filename.c_str());
	for (int j=0; j<ny; ++j){
		for (int i=0; i<nx; ++i){
			fout << elev[j*nx+i] << "\t";
		}
		fout << "\n";
	}
	fout.close();
}





void PointCloud::read_las(string file){
	ifstream ifs;
	ifs.open(file.c_str(),ios::in | ios::binary);
	if (!ifs) cout << "Error Opening file: " << file << endl;

	liblas::ReaderFactory f;
	liblas::Reader reader = f.CreateWithStream(ifs);
	liblas::Header const& header = reader.GetHeader();
	cout << "Compressed: " << ((header.Compressed() == true) ? "true\n":"false\n");
	cout << "Signature: " << header.GetFileSignature() << '\n';
	cout << "Points count: " << header.GetPointRecordsCount() << '\n';
	nverts = header.GetPointRecordsCount();

	points.reserve(3*nverts);
	
	while (reader.ReadNextPoint()){ //count<10
		liblas::Point const& p = reader.GetPoint();
		points.push_back(p.GetX());
		points.push_back(p.GetY());
		points.push_back(p.GetZ());

		//cout << p.GetX() << ", " << p.GetY() << ", " << p.GetZ() << "\n";
	}

}


void PointCloud::createDEM(float dx, float dy){
	float xmin = min_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), 
							 [](const fl::vec3& p1, const fl::vec3& p2){return (p1.x < p2.x);} )->x;
	float xmax = max_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), 
							 [](const fl::vec3& p1, const fl::vec3& p2){return (p1.x < p2.x);} )->x;
	float ymin = min_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), 
							 [](const fl::vec3& p1, const fl::vec3& p2){return (p1.y < p2.y);} )->y;
	float ymax = max_element((fl::vec3*)points.data(), (fl::vec3*)(points.data()+3*nverts), 
							 [](const fl::vec3& p1, const fl::vec3& p2){return (p1.y < p2.y);} )->y;
	cout << "x range: " << xmin << " " << xmax << endl;	  
	cout << "y range: " << ymin << " " << ymax << endl;		
	
	xmin = floor(xmin/dx)*dx;
	xmax = ceil(xmax/dx)*dx;
	ymin = floor(ymin/dy)*dy;
	ymax = ceil(ymax/dy)*dy;

	dem.initFromBounds(xmin,xmax,ymin,ymax,dx,dy);
	//dem.init(nx,ny,dx,dy,xmin,ymin);
	
	fill(dem.elev.begin(), dem.elev.end(), 1e20);	// initialize all elevations with 1e20
	for (int k=0; k<nverts; ++k){
		int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
		int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
		dem.elev[iy*dem.nx+ix] = min(points[3*k+2], dem.elev[iy*dem.nx+ix]);
	}
	
	// TODO: Filter DEM
}

void PointCloud::subtractDEM(){
	for (int k=0; k<nverts; ++k){
		int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
		int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
		points[3*k+2] -= dem.elev[iy*dem.nx+ix];	// TODO: subract interpolated value
	}
}	

void PointCloud::deleteGround(float ht){
	for (int k=0; k<nverts; ++k){
		int ix = floor((points[3*k+0]- dem.xmin)/ dem.dx);
		int iy = floor((points[3*k+1]- dem.ymin)/ dem.dy);
		if (points[3*k+2] < ht) points[3*k+0] = points[3*k+1] = points[3*k+2] = 0; 
	}
}

void PointCloud::generateRandomClusters(int N, int nc, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float R){
	nverts = N;
	points.resize(3*nverts);
	
	// create centroids
	vector <float> centroids;
	centroids.reserve(nc);
	for (int i=0; i<nc; ++i){
		centroids.push_back(runif(xmin, xmax));
		centroids.push_back(runif(ymin, ymax));
		centroids.push_back(runif(zmin, zmax));
	}

	for (int i=0; i<N; ++i){
		int g = rand() % nc;
		float dx = runif(-R, R);
		float dy = runif(-R, R);
		float dz = runif(-R, R);
	
		points[3*i+0] = centroids[3*g+0] + dx;
		points[3*i+1] = centroids[3*g+1] + dy;
		points[3*i+2] = centroids[3*g+2] + dz;
	}
}





