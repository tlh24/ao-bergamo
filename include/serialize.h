// class for serializing (matlab format) to a file
// & read (binary format) from a mmaped file.
#ifndef __SERIALIZE_H__
#define __SERIALIZE_H__
#include <string>
#include <vector>
#include <array>

using namespace std;
class Serialize
{
public:
	string	m_name;
	int		m_lastBackup; //index+1 of last record saved to backup.
	Serialize() {}
	virtual ~Serialize() {
		m_lastBackup = 0;
	}
	void perr(const char *method) {
		fprintf(stderr, "\"%s\":%s must be implemented in derived classes.\n",
		        m_name.c_str(), method);
	}
	virtual bool 	store() {
		perr("store");
		return false;
	}
	virtual void	clear() {
		perr("clear");
	}
	virtual int 	nstored() {
		return 0;       //number of timeslices.
	}
	virtual string storeName(int ) {
		perr("storeName");
		return string("none");
	}
	virtual int 	getStoreClass(int ) {
		perr("getStoreClass");
		return 0;
	}
	virtual void	getStoreDims(int , size_t *dims) {
		perr("getStoreDims");
		dims[0] = 0;
		dims[1] = 0;
	}
	virtual void	*getStore(int , int) {
		perr("getStore");
		return NULL;
	}
	virtual int 	numStores() {
		perr("numStores");
		return 0;
	}
	virtual double 	*mmapRead(double * ) {
		perr("mmapRead");
		return NULL;
	}
	//drawing routines -- opengl -- not all need implement.
	virtual void	draw(int, float) {}
	virtual void	move(long double) {}
	// reads/writes parameters from a mmaped file (address).
	// all mmap variables are doubles, for convenience.
	virtual string getMmapInfo() {
		std::stringstream oss;
		size_t dims[2];
		for (int indx = 0; indx < numStores(); indx++) {
			getStoreDims(indx, dims);
			string stor = storeName(indx);
			oss << "\t'double' [" << dims[0] << " " << dims[1] << "] '" << stor << "';...\n";
		}
		return oss.str();
	}
	virtual string getStructInfo() {
		std::stringstream oss;
		size_t dims[2];
		for (int indx = 0; indx < numStores(); indx++) {
			getStoreDims(indx, dims);
			string stor = storeName(indx);
			oss << "b5." << stor << " = zeros(" << dims[0] << "," << dims[1] << ");\n";
		}
		return oss.str();
	}
};


// this class is for recording arbitrary numbers, in the form of a vector,
// from matlab.  Can be used for e.g. trial#, trial type, your mommas number ...
// ... pager number, of course.
template <class T>
class VectorSerialize : public Serialize
{
public:
	int	m_size;
	int	m_type;
	double				m_time;
	vector<T> 			m_stor;
	vector<double>		v_time;
	vector<vector<T> > v_stor;
	T		*m_bs;

	VectorSerialize(int size, int matiotype) : Serialize() {
		m_size = size;
		m_type = matiotype;
		for (int i=0; i<size; i++) {
			m_stor.push_back((T)0);
		}
		m_bs = NULL;
	}
	~VectorSerialize() {
		clear();
		free(m_bs);
	}
	virtual bool store() {
		bool same = true; //delta compression.
		int n = nstored();
		if (n > 0) {
			for (int i=0; i<m_size; i++)
				same &= (m_stor[i] == v_stor[n-1][i]);
		} else same = false;
		if (!same) {
			m_time = gettime();
			v_time.push_back(m_time);
			v_stor.push_back(m_stor);
		}
		return !same;
	}
	virtual void clear() {
		v_time.clear();
		v_stor.clear();
	}
	virtual int nstored() {
		return v_stor.size();
	}
	virtual string storeName(int indx) {
		switch (indx) {
		case 0:
			return m_name + string("time_o");
		case 1:
			return m_name + string("v");
		}
		return string("none");
	}
	virtual int getStoreClass(int indx) {
		switch (indx) {
		case 0:
			return MAT_C_DOUBLE;
		case 1:
			return m_type;
		}
		return 0;
	}
	virtual void getStoreDims(int indx, size_t *dims) {
		switch (indx) {
		case 0:
			dims[0] = 1;
			dims[1] = 1;
			return;
		case 1:
			dims[0] = m_size;
			dims[1] = 1;
			return;
		}
	}
	virtual void *getStore(int indx, int k) {
		//coalesce the memory -- <vector<vector>> is non-continuous in memory.
		if (indx == 0) {
			return (void *)&(v_time[k]);
		} else if (indx == 1) {
			if (m_bs) free(m_bs);
			int n = nstored(); //atomic -- if we're not careful, may change during read!
			m_bs = (T *)malloc(sizeof(T)*(n-k)*m_size);
			for (int i=0; i<n-k; i++) {
				for (int j=0; j<m_size; j++) {
					m_bs[j + i*m_size] = v_stor[i+k][j];
				}
			}
			return (void *)(m_bs);
		} else return NULL;
	}
	virtual int numStores() {
		return 2;
	}
	virtual double *mmapRead(double *d) {
		*d++ = m_time; //when the vector was updated.
		for (int i=0; i<m_size; i++) {
			m_stor[i] = (T)(*d++); //default input.
		}
		return d;
	}
};

// same as VectorSerialize, except two variables -- position and velocity e.g.
// yea, copypasta. sorry.
template <class T>
class VectorSerialize2 : public Serialize
{
public:
	int	m_size;
	int	m_type;
	double				m_time;
	vector<T> 			m_stor;
	vector<T> 			m_stor2;
	vector<double>		v_time;
	vector<vector<T> > v_stor;
	vector<vector<T> > v_stor2;
	T		*m_bs;

	VectorSerialize2(int size, int matiotype) : Serialize() {
		m_size = size;
		m_type = matiotype;
		for (int i=0; i<size; i++) {
			m_stor.push_back((T)0);
			m_stor2.push_back((T)0);
		}
		m_bs = NULL;
	}
	~VectorSerialize2() {
		clear();
		free(m_bs);
	}
	virtual bool store() {
		bool same = true; //delta compression.
		int n = nstored();
		if (n > 0) {
			for (int i=0; i<m_size; i++) {
				same &= (m_stor[i] == v_stor[n-1][i]);
				same &= (m_stor2[i] == v_stor2[n-1][i]);
			}
		} else same = false;
		if (!same) {
			m_time = gettime();
			v_time.push_back(m_time);
			v_stor.push_back(m_stor);
			v_stor2.push_back(m_stor2);
		}
		return !same;
	}
	virtual void clear() {
		v_time.clear();
		v_stor.clear();
		v_stor2.clear();
	}
	virtual int nstored() {
		return v_stor.size();
	}
	virtual string storeName(int indx) {
		switch (indx) {
		case 0:
			return m_name + string("time_o");
		case 1:
			return m_name + string("x");
		case 2:
			return m_name + string("y");
		}
		return string("none");
	}
	virtual int getStoreClass(int indx) {
		switch (indx) {
		case 0:
			return MAT_C_DOUBLE;
		case 1:
		case 2:
			return m_type;
		}
		return 0;
	}
	virtual void getStoreDims(int indx, size_t *dims) {
		switch (indx) {
		case 0:
			dims[0] = 1;
			dims[1] = 1;
			return;
		case 1:
		case 2:
			dims[0] = m_size;
			dims[1] = 1;
			return;
		}
	}
	virtual void *getStore(int indx, int k) {
		//coalesce the memory -- <vector<vector>> is non-continuous in memory.
		//k is the starting offset.  caller must free the memory.
		if (indx == 0) {
			return (void *)&(v_time[k]);
		} else if (indx == 1 || indx == 2) {
			if (m_bs) free(m_bs);
			int n = nstored(); //atomic -- if we're not careful, may change during read!
			m_bs = (T *)malloc(sizeof(T)*(n-k)*m_size);
			for (int i=0; i<n-k; i++) {
				if (indx == 1) {
					for (int j=0; j<m_size; j++)
						m_bs[j + i*m_size] = v_stor[i+k][j];
				} else {
					for (int j=0; j<m_size; j++)
						m_bs[j + i*m_size] = v_stor2[i+k][j];
				}
			}
			return (void *)(m_bs);
		} else return NULL;
	}
	virtual int numStores() {
		return 3;
	}
	virtual double *mmapRead(double *d) {
		*d++ = m_time; //when the vector was updated.
		for (int i=0; i<m_size; i++) {
			m_stor[i] = (T)(*d++);
		}
		for (int i=0; i<m_size; i++) {
			m_stor2[i] = (T)(*d++);
		}
		return d;
	}
};

//convenience class for saving 4x4 calibration matrix (affine, quadratic).
class Matrix44Serialize : public Serialize
{
public:
	array<double,16>			m_cmp;
	array<float,16>				m_x;
	vector<array<float,16> >	v_x;

	Matrix44Serialize(string name) : Serialize() {
		m_name = name;
		for (int i=0; i<16; i++) {
			m_x[i] = 0; // Matlab ordering (column major -- same in openGL).
			m_cmp[i] = 0;
		}
		for (int i=0; i<4; i++) {
			m_x[i+i*4] = 1.f;
			m_cmp[i+i*4] = 1.f;
		}
	}
	~Matrix44Serialize() {
		clear();
	}
	virtual void clear() {
		v_x.clear();
	}
	virtual bool store() {
		return false; /* in mmap. */
	}
	virtual int nstored() {
		return v_x.size();
	}
	virtual string storeName(int ) {
		return m_name + string("_m44");
	}
	virtual int getStoreClass(int ) {
		return MAT_C_SINGLE;
	}
	virtual void getStoreDims(int, size_t *dims) {
		dims[0] = 4;
		dims[1] = 4;
		return;
	}
	virtual void *getStore(int , int k) {
		return (void *)&(v_x[k]);
	}
	virtual int numStores() {
		return 1;
	}
	virtual double *mmapRead(double *d) {
		//don't store the time here, as an incentive to not change it during the exp!
		bool sames = true;
		for (int i=0; i<16; i++)
			sames &= (d[i] == m_cmp[i]);
		if (!sames) {
			for (int i=0; i<16; i++) {
				m_x[i] = (float)d[i];
				m_cmp[i] = d[i];
			}
			v_x.push_back(m_x);
		}
		d += 16;
		return d;
	}
	float *data() {
		return m_x.data();
	};
};

// class MouseSerialize : public VectorSerialize<float>
// {
// public:
// 	//timing via VectorSerialize.
// 	MouseSerialize() : VectorSerialize(3, MAT_C_SINGLE) {
// 		m_name = "mouse_";
// 	}
// 	~MouseSerialize() {}
// 	virtual bool store() {
// 		m_stor[0] = g_mousePos[0];
// 		m_stor[1] = g_mousePos[1];
// 		m_stor[2] = 0.f;
// 		return VectorSerialize::store();
// 	}
// 	virtual string storeName(int ) {
// 		return m_name + string("sensors_o");        //output.
// 	}
// 	virtual double *mmapRead(double *d) {
// 		*d++ = g_mousePos[0];
// 		*d++ = g_mousePos[1];
// 		*d++ = 0.f;
// 		return d;
// 	}
// };


void writeMatlab(vector<Serialize *> tosave, char *filename, bool backup);
size_t matlabFileSize(vector<Serialize *> tosave);
size_t mmapFileSize(vector<Serialize *> tosave);
bool matlabHasNewData(vector<Serialize *> tosave);
#endif

