#include "matrix.h"
#include <assert.h>
#include <algorithm>
#include "myfunctions.h"
#include "testing_parameters.h"
#include "boundaries.h"
#include "boundary_conditions.h"
#include <iostream>

using namespace myvector;
using namespace std;
using namespace testingparameters;
using namespace boundary_conditions;
using namespace partition;
using namespace boundaries;
using namespace logger;

namespace matrix
{	
	Matrix::Matrix(){};

	Matrix::Matrix(int size1, int size2)
	{
		initialize(size1, size2);
	}

	void Matrix::initialize(int size1, int size2)
	{		
		n = size1; size = size2;

		initialize_vector(ggl, size);
		initialize_vector(ggu, size);
		initialize_vector(di, n);
		initialize_vector(ig, n + 1);
		initialize_vector(jg, size);
		yl.initialize(n);
		yu.initialize(n);

		if (Testing_parameters::use_LU)
		{
			initialize_vector(LU_ggl, size);
			initialize_vector(LU_ggu, size);
			initialize_vector(LU_di, n);
		}
	};

	void Matrix::reinitialize()
	{
		memset(&ggl[0], 0, size * sizeof(double)); //обнуляем
		memset(&ggu[0], 0, size * sizeof(double)); //обнуляем
		memset(&di[0], 0, n * sizeof(double)); //обнуляем

		if (Testing_parameters::use_LU)
		{
			memset(&LU_ggl[0], 0, size * sizeof(double));
			memset(&LU_ggu[0], 0, size * sizeof(double));
			memset(&LU_di[0], 0, n * sizeof(double));
		}
	};

	MyVector Matrix::operator*(MyVector& a) 
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(a.ar.size());

		assert(a.ar.size() == n);
		for(i = 0; i < n; i++)
		{
			kol = ig[i + 1] - ig[i];//количество ненулевых элементов строки (столбца) от первого
									//ненулевого элемента до диагонального элемента (не включая его)
			iend = ig[i + 1];
			k = ig[i]; // адрес первого занятого элемента строки (столбца) 

			new_vector[i] = di[i] * a[i];//от главной диагонали

			for(; k < iend; k++)//проходим по всем элементам i строки (столбца)
			{
				j = jg[k];
				new_vector[i] += ggl[k] * a[j];//от нижнего треугольника
				new_vector[j] += ggu[k] * a[i];//от верхнего треугольника
			}
		}

		return new_vector;
	}

	MyVector Matrix::operator/(MyVector& a) 
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(a.ar.size());

		assert(a.ar.size() == n);
		for(i = 0; i < n; i++)
		{
			kol = ig[i + 1] - ig[i];//количество ненулевых элементов строки (столбца) от первого
									//ненулевого элемента до диагонального элемента (не включая его)
			iend = ig[i + 1];
			k = ig[i]; // адрес первого занятого элемента строки (столбца) 

			new_vector[i] = di[i] * a[i];//от главной диагонали

			for(; k < iend; k++)//проходим по всем элементам i строки (столбца)
			{
				j = jg[k];
				new_vector[i] += ggu[k] * a[j];//от нижнего треугольника
				new_vector[j] += ggl[k] * a[i];//от верхнего треугольника
			}
		}

		return new_vector;
	}

	Matrix::~Matrix(){};

	MyVector Matrix::Uv(MyVector& v)
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(v.ar.size());

		assert(v.ar.size() == n);
		for(i = 0; i < n; i++)
		{
			kol = ig[i+1] - ig[i];//количество ненулевых элементов столбца от первого
									//ненулевого элемента до диагонального элемента (не включая его)
			iend = ig[i+1];
			k = ig[i]; // адрес первого занятого элемента столбца

			new_vector[i] = v[i];//от главной диагонали (у U на диагонали 1)

			for(; k < iend; k++)//проходим по всем элементам i столбца
			{
				j = jg[k];
				new_vector[j] += ggu[k] * v[i];//от верхнего треугольника
			}
		}

		return new_vector;
	}

	void Matrix::add_element(int i, int j, double element)
	{
		int id;
		bool flag;

		if(i == j)
			di[i] += element;
		else
		{
			if(i < j)
			{	
				flag = false;
				for(id = ig[j]; !flag && id <= ig[j + 1] - 1; id++)
					if(jg[id] == i) flag = true;
				 if(flag) ggu[id - 1] += element;
			}
			else
			{
				flag = false;
				for(id = ig[i]; !flag && id <= ig[i + 1] - 1; id++)
					if(jg[id] == j) flag = true;
				if(flag) ggl[id - 1] += element;
			}
		}
	}

	void Matrix::put_element(int i, int j, double element)
	{
		int id;
		bool flag;

		if(i == j)
			di[i] = element;
		else
		{
			if(i < j)
			{	
				flag = false;
				for(id = ig[j]; !flag && id <= ig[j + 1] - 1; id++)
					if(jg[id] == i) flag = true;
				 if(flag) ggu[id - 1] = element;
			}
			else
			{
				flag = false;
				for(id = ig[i]; !flag && id <= ig[i + 1] - 1; id++)
					if(jg[id] == j) flag = true;
				if(flag) ggl[id - 1] = element;
			}
		}
	}

	void Matrix::create_portret(Partition& p, Logger& log)
	{
		log.send_message_create_portret();

		vector <int> unzero_elements_list;
		//vector <int> *lists;
		vector<vector <int>> lists; lists.resize(n);
		int unzero_elements_lists_size;
		int current_number;

		//lists = new vector <int>[n];
		int count_elements = p.elements.size();

		for (int i = 0; i < count_elements; i++)
		{
			//собираем ненулевые для к.э. глобальные номера узлов
			unzero_elements_lists_size = p.create_unzero_elements_list(i, unzero_elements_list);
			//т.к. первые 9 штук в списке - номера узлов элемента,то 
			//идём по ним и выбираем для каждого номера, меньшие его,
			//т.е. те, которые будут располагаться левее соответствующей
			//диагонали, потому что портрет строим по строкам
			//а затем кладём в соответствующий список
			int nf = p.elements[i].ndof;
			for (int j = 0; j < nf; j++)
			{
				current_number = unzero_elements_list[j];
				for (int k = 0; k < unzero_elements_lists_size; k++)
					if (unzero_elements_list[k] < current_number)
					{ 
						lists[current_number].push_back(unzero_elements_list[k]);
					}
				sort(lists[current_number].begin(), lists[current_number].end());
			}
			unzero_elements_list.clear();
		}

		ig[0] = 0;

		for (int i = 0; i < n; i++)
		{
			if (!lists[i].empty())
				ig[i + 1] = ig[i] + lists[i].size();
			else ig[i + 1] = ig[i];
		}

		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (!lists[i].empty())
			{
				for (unsigned int j = 0; j < lists[i].size(); j++)
				{
					jg[k] = lists[i][j];
					k++;
				}
				lists[i].clear();
			}
		}

		//delete[] lists;
	}
}