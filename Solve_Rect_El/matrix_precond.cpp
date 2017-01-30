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
	void Matrix::LU()
	{
		int i;
		int i0, j0;
		int iend;
		int num;
		int ki, kj;
		double suml, sumu, sumdg;
		int size2 = size;

		for (i = 0; i < size2; i++)
		{
			LU_ggu[i] = ggu[i];
			LU_ggl[i] = ggl[i];
		}

		for (i = 0; i < n; i++)
			LU_di[i] = di[i];

		for (i = 0; i < n; i++)
		{
			i0 = ig[i];
			iend = ig[i + 1];
			for (num = i0, sumdg = 0; num < iend; num++)
			{
				j0 = ig[jg[num]]; //� ����������� �� ������ ��������� �������,����� ������� l,������ �������  ���� ��������� �� � u 
				int jend = ig[jg[num] + 1];
				ki = i0;
				kj = j0;
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++) //��� num ����������� ��� ���������� ��������
				{
					for (int m = kj; m < jend; m++)
						if (jg[ki] == jg[m]) //���� ��������������� ��������� �������� ��� ���������
						{
							suml += LU_ggl[ki] * LU_ggu[m];
							sumu += LU_ggl[m] * LU_ggu[ki];//��� ������������� �������� �� U
						}
				}
				LU_ggl[num] -= suml;
				LU_ggu[num] = (LU_ggu[num] - sumu) / LU_di[jg[num]];
				sumdg += LU_ggl[num] * LU_ggu[num];//���������� ������������ ��������	
			}
			LU_di[i] -= sumdg;
		}
		//int a, b, j;
		//LU_di[0] = di[0];
		//for (i = 0; i<n; i++)
		//{
		//	LU_di[i] = di[i];

		//	for (j = ig[i]; j<ig[i + 1]; j++)

		//	{
		//		LU_ggl[j] = ggl[j];
		//		LU_ggu[j] = ggu[j];

		//		for (a = ig[i]; a<j; a++)
		//		{
		//			for (b = ig[jg[j]] - 1; b <= ig[jg[j]+1] - 1; b++)
		//			{
		//				if (jg[a] == jg[b])
		//				{
		//					LU_ggl[j] -= LU_ggl[a] * LU_ggu[b];
		//					LU_ggu[j] -= LU_ggu[a] * LU_ggl[b];
		//				}
		//			}
		//		}
		//		LU_ggu[j] /= LU_di[jg[j]];
		//		LU_di[i] -= LU_ggu[j] * LU_ggl[j];
		//	}
		//}
	}
	void Matrix::LYF(MyVector& b)
	{
		int i, k;
		int i0;//����� ������ ������
		int iend;//����� ����� ������
		double sum;

		yl.make_zero();

		if (Testing_parameters::use_LU)
		{
			for (i = 0; i < n; i++)
			{
				i0 = ig[i]; iend = ig[i + 1];

				for (k = i0, sum = 0; k < iend; k++)
					sum += yl[jg[k]] * LU_ggl[k];

				yl[i] = (b[i] - sum) / LU_di[i];
			}
		}
		else
			yl = b;

		//int i, k;

		//yl.make_zero();

		//if (Testing_parameters::use_LU)
		//{
		//	for (i = 0; i < n; i++)
		//	{
		//		yl[i] = b[i];

		//		for (k = ig[i]; k < ig[i + 1]; k++)
		//			yl[i] -= yl[jg[k]] * LU_ggl[k];

		//		yl[i] /= LU_di[i];
		//	}
		//}
		//else
		//	yl = b;
	}
	void Matrix::LYFt(MyVector& b)
	{
		int i, k;
		int i0;//����� ������ ������
		int iend;//����� ����� ������
		double sum;

		yl.make_zero();
		if (Testing_parameters::use_LU)
		{
			MyVector bb(n);
			bb = b;
			for (i = n - 1; i >= 0; i--)
			{
				i0 = ig[i]; iend = ig[i + 1];
				yl[i] = bb[i] / LU_di[i];
				for (k = i0, sum = 0; k < iend; k++)
					bb[jg[k]] -= yl[i] * LU_ggl[k];
			}
		}
		else
			yl = b;

		//int i, k;
		//MyVector bb(n); bb = b;

		//yl.make_zero();
		//if (Testing_parameters::use_LU)
		//{
		//	MyVector bb(n);
		//	bb = b;
		//	for (i = n - 1; i >= 0; i--)
		//	{
		//		yl[i] = bb[i] / LU_di[i];
		//		for (k = ig[i]; k < ig[i + 1]; k++)
		//			bb[jg[k]] -= yl[i] * LU_ggl[k];
		//	}
		//}
		//else
		//	yl = b;
	}
	void Matrix::UXY(MyVector& b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if (Testing_parameters::use_LU)
		{
			for (i = n - 1; i >= 0; i--)//������ �� �������� � �����
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for (k = iend - 1; k >= i0; k--)//��� �� ������� � �����
					yu[jg[k]] -= yu[i] * LU_ggu[k];
			}
		}
		else
			yu = b;

		//int i, k;

		//MyVector bb(n); bb = b;
		//yu.make_zero();
		//if (Testing_parameters::use_LU)
		//{
		//	for (i = n - 1; i >= 0; i--)//������ �� �������� � �����
		//	{
		//		yu[i] = bb[i];

		//		for (k = ig[i]; k < ig[i + 1]; k++)//��� �� ������� � �����
		//			bb[jg[k]] -= yu[i] * LU_ggu[k];
		//	}
		//}
		//else
		//	yu = b;
	}
	void Matrix::UXYt(MyVector& b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if (Testing_parameters::use_LU)
		{
			for (i = n - 1; i >= 0; i--)//������ �� �������� � �����
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for (k = iend - 1; k >= i0; k--)//��� �� ������� � �����
					yu[i] -= yu[jg[k]] * LU_ggu[k];
			}
		}
		else
			yu = b;

		//int i, k;

		//yu.make_zero();
		//if (Testing_parameters::use_LU)
		//{
		//	for (i = 0; i < n; i++)//������ �� �������� � �����
		//	{
		//		yu[i] = b[i];

		//		for (k = ig[i]; k < ig[i + 1]; k++)//��� �� ������� � �����
		//			yu[i] -= yu[jg[k]] * LU_ggu[k];
		//	}
		//}
		//else
		//	yu = b;
	}

	MyVector Matrix::Uv_(MyVector& v)
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(v.ar.size());

		assert(v.ar.size() == n);
		if (Testing_parameters::use_LU)
		{
			for (i = 0; i < n; i++)
			{
				kol = ig[i + 1] - ig[i];//���������� ��������� ��������� ������� �� �������
										//���������� �������� �� ������������� �������� (�� ������� ���)
				iend = ig[i + 1];
				k = ig[i]; // ����� ������� �������� �������� �������

				new_vector[i] = v[i];//�� ������� ��������� (� U �� ��������� 1)

				for (; k < iend; k++)//�������� �� ���� ��������� i �������
				{
					j = jg[k];
					new_vector[j] += LU_ggu[k] * v[i];//�� �������� ������������
				}
			}
		}
		else
			new_vector = v;
		return new_vector;
	}

}