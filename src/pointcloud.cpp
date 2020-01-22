#include "../headers/pointcloud.h"
#include <cmath>
#include <liblas/liblas.hpp>
#include "../utils/simple_math.h"

using namespace std;


bool compare_z(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.z < p2.z);
}

bool compare_y(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.y < p2.y);
}

bool compare_x(const fl::vec3& p1, const fl::vec3& p2){
	return (p1.x < p2.x);
}

void sort_by_y(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_y); 
}
	
void sort_by_z(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_z); 
}

void sort_by_x(float* begin, float* end){
	sort((fl::vec3*)begin, (fl::vec3*)end, compare_x); 
}




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

void PointCloud::write_las(string file){

	std::ofstream ofs;
	ofs.open(file.c_str(), ios::out | ios::binary);

	liblas::Header header;
//	header.SetDataFormatId(liblas::ePointFormat1); // Time only

//	// Set coordinate system using GDAL support
//	liblas::SpatialReference srs;
//	srs.SetFromUserInput("EPSG:4326");

//	header.SetSRS(srs);
//	header.setPointRecordsCount(nverts);

	liblas::Writer writer(ofs, header);
	// here the header has been serialized to disk into the *file.las*

	liblas::Point point(&header);
	for (int i=0; i<nverts; ++i){
		point.SetCoordinates(points[3*i+0], points[3*i+1], points[3*i+2]);
		writer.WritePoint(point);
	}
	
	ofs.close();
}


void PointCloud::write_ascii(string file){

	ofstream fout(file.c_str());
	for (int i=0; i< nverts; ++i){
		fout << points[3*i+0] << "\t" << points[3*i+1] << "\t" << points[3*i+2] << "\n";
	}
	fout.close();
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

void PointCloud::subtractDEM_bil(){
	for (int k=0; k<nverts; ++k){
		float x = points[3*k+0], y = points[3*k+1];
		int ix = floor((x - dem.xmin)/ dem.dx);
		int iy = floor((y - dem.ymin)/ dem.dy);

		float x1 = dem.xmin + ix*dem.dx;
		float x2 = dem.xmin + (ix+1)*dem.dx;
		float y1 = dem.ymin + iy*dem.dy;
		float y2 = dem.ymin + (iy+1)*dem.dy;

		float f0 = dem.elev[iy*dem.nx+ix];
		float f1 = dem.elev[(iy+1)*dem.nx+ix];
		float f2 = dem.elev[iy*dem.nx+ix+1];
		float f3 = dem.elev[(iy+1)*dem.nx+ix+1];

		if (f1 > 1e19) f1 = f0;
		if (f2 > 1e19) f2 = f0;
		if (f3 > 1e19) f3 = f0;

		float f = (f0 * ((x2-x )*(y2-y ))   // f11
				 + f2 * ((x -x1)*(y2-y ))   // f21
				 + f1 * ((x2-x )*(y -y1))   // f12
				 + f3 * ((x -x1)*(y -y1)))  // f22
				 / ((x2-x1)*(y2-y1));		// x2=x1+dx, so no div-by-zero!

//		if (k % 1000 == 0) cout << "elev corners: " << f0 << " " << f1 << " " <<  f2 << " " << f3 << "|" << f << endl;
		
		points[3*k+2] -= f;	// TODO: subract interpolated value

//		if (k % 1000 == 0) cout << points[3*k+0] << " " << points[3*k+1] << " " << points[3*k+2] << endl;
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



void PointCloud::calcGridParams(float Rg, Grid & par){
	float3 * pos = (float3*)points.data();
	par.cellSize.x = par.cellSize.y = par.cellSize.z = Rg;
	
	float xmin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_x)->x;
	float xmax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_x)->x;
	float ymin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_y)->y;
	float ymax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_y)->y;
	float zmin = min_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_z)->z;
	float zmax = max_element((fl::vec3*)pos, (fl::vec3*)(pos+nverts), compare_z)->z;

	cout << "Points:\n";
	cout << "x range: " << xmin << " " << xmax << endl;	  
	cout << "y range: " << ymin << " " << ymax << endl;		
	cout << "z range: " << zmin << " " << zmax << endl;		
	
	xmin = floor(xmin/Rg)*Rg;
	xmax = ceil(xmax/Rg)*Rg;
	ymin = floor(ymin/Rg)*Rg;
	ymax = ceil(ymax/Rg)*Rg;
	zmin = floor(zmin/Rg)*Rg;
	zmax = ceil(zmax/Rg)*Rg;

	float3 origin;
	origin.x = xmin; origin.y = ymin; origin.z = zmin;
	par.origin = origin;

	par.gridSize.x  =  (xmax- xmin)/ Rg ; 
	par.gridSize.y  =  (ymax- ymin)/ Rg ; 
	par.gridSize.z  =  (zmax- zmin)/ Rg ; 

	cout << "Grid (spacing = " << Rg << ")\n";
	cout << "x range: " << xmin << " " << xmax << endl;	  
	cout << "y range: " << ymin << " " << ymax << endl;		
	cout << "z range: " << zmin << " " << zmax << endl;		

	cout << "gridSize = (" << par.gridSize.x << ", " << par.gridSize.y << ", " << par.gridSize.z << ")" << endl;
}


float PointCloud::distance(int p, int q){
	float dx = points[3*p+0] - points[3*q+0];
	float dy = points[3*p+1] - points[3*q+1];
	float dz = points[3*p+2] - points[3*q+2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}



