#include "suAgent.h"
#include "suStructureEvolution.h"
#include "suVolume.h"
#include "suMorton.h"
#include <deque>
#include "config.h"
#include "suWriteSTL.h"


#include <Eigen/Core>

int suStructrueOptimizer::init_morton_coding()
{
	nodeArr_.resize(1 << (3 * nLevel_));
	std::vector<SU::OctNode*> &nodePointers = pVolume_->leafBoundaryNodes_;

	SU::Point minP(pVolume_->bbMin_[0], pVolume_->bbMin_[1], pVolume_->bbMin_[2]);

	for (int i = 0; i < (int)nodePointers.size(); i++)
	{
		if (nodePointers[i]->label_ == SU::INTERIOR_CELL||nodePointers[i]->label_==SU::BOUNDARY_CELL||nodePointers[i]->label_==SU::BOUNDARY_CELL_SPECIAL)
		{
			nodePointers[i]->out = true;
			SU::Point _cp = nodePointers[i]->center() - minP;
			SU::Point _vsize = nodePointers[i]->size();
			nodePointers[i]->morton = suMorton::encode(_cp.x / _vsize.x, _cp.y / _vsize.y, _cp.z / _vsize.z, nLevel_);
			nodeArr_[nodePointers[i]->morton] = nodePointers[i];
		}
	}

	nodePointers = pVolume_->leafInternalNodes_;
	for (int i = 0; i < (int)nodePointers.size(); i++)
	{
		if (nodePointers[i]->label_ == SU::INTERIOR_CELL || nodePointers[i]->label_ == SU::BOUNDARY_CELL || nodePointers[i]->label_ == SU::BOUNDARY_CELL_SPECIAL)
		{
			nodePointers[i]->out = true;
			SU::Point _cp = nodePointers[i]->center() - minP;
			SU::Point _vsize = nodePointers[i]->size();
			nodePointers[i]->morton = suMorton::encode(_cp.x / _vsize.x, _cp.y / _vsize.y, _cp.z / _vsize.z, nLevel_);
			nodeArr_[nodePointers[i]->morton] = nodePointers[i];
		}
	}
	
	/*int count = 0;
	auto countIt = nodeArr_.begin();
	for (; countIt != nodeArr_.end(); countIt++) {
		if ((*countIt) == NULL)
			if ((*countIt)->label_ == SU::INTERIOR_CELL)
				count++;
	}
	std::cout << count; system("pause");*/
	return (int)nodeArr_.size();
}
int suStructrueOptimizer::init_position_field()
{
	std::deque<SU::OctNode*> floodEdges;
	int maxValue = 0;

	///generate nodes list(extending all nodes to leaf nodes with min size)	
	//std::vector<SU::OctNode*> &BorderNodeArr = pVolume_->leafBoundaryNodes_;
	//nodeArr_ = pVolume_->leafBoundaryNodes_;   //copy node pointer


	//init boundary
	for (int i = 0; i < (int)nodeArr_.size(); i++)
	{
		if (nodeArr_[i] == NULL) {
			nodeArr_[i] = new SU::OctNode;
			nodeArr_[i]->label_ = SU::EXTERIOR_CELL;
			nodeArr_[i]->morton = i;
		}
		if (nodeArr_[i]->label_ == SU::BOUNDARY_CELL || nodeArr_[i]->label_ == SU::BOUNDARY_CELL_SPECIAL)
		{
			nodeArr_[i]->location = 1;
			floodEdges.push_back(nodeArr_[i]);
		}

	}


	//floodfill from boundary to internal
	int ii = 0;
	while (!(floodEdges.empty()))
	{
		std::deque<SU::OctNode*>::iterator it = floodEdges.begin();

		std::vector<SU::OctNode*> neighbors;  //use octree to find neighbor, because the searching is from boundary nodes
		SU::OctNode *pCurNode = *it;

		pVolume_->get6NeighborNodes(*it, neighbors);
		for (int i = 0; i < (int)neighbors.size(); i++)
		{
			/*int mcode = neighbors[i];
			if (mcode >= (int)nodeArr_.size() || mcode < 0) continue;
			SU::OctNode *pNode = nodeArr_[mcode];*/
			SU::OctNode *pNode = neighbors[i];
			if (pNode)
				if (pNode->location == 0 &&
					pNode->label_ == SU::INTERIOR_CELL)
				{
					maxValue = pCurNode->location + 1;
					pNode->location = maxValue;
					floodEdges.push_back(pNode);
				}
		}

		floodEdges.pop_front();
	}

	return maxValue;

}

int suStructrueOptimizer::init(SU::suVolume *pVolume)
{
	pVolume_ = pVolume;
	init_morton_coding();
	init_position_field();

	return 0;
}



//�ڱ����Ƴ����£�����ʴ����������Ҫ����Ϊ����Ϊ��ʵ���ȶ������ƣ������͡��ǲ����ٵġ�

int suStructrueOptimizer::evolve()
{
	SU::suAgent<SU::OctNode> agent;

	std::vector<SU::suAgent<SU::OctNode> > agents;

	for (unsigned int i = 0; i < nodeArr_.size(); i++)
	{
		if (nodeArr_[i]->label_ == SU::EXTERIOR_CELL)
			continue;
		std::vector<int> neighbors;
		std::vector<SU::OctNode *> env;
		suMorton::get_26neighbors(neighbors, nodeArr_[i]->morton, nodeArr_[i]->level_);
		for (int j = 0; j < neighbors.size(); j++)
		{
			env.push_back(nodeArr_[neighbors[j] ] );
		}
		agent.bind(nodeArr_[i], env);
		agents.push_back(agent);
	}	

	int nIteration = 0;

	/*while (true)
	{*/
	for (unsigned int i = 0; i < agents.size(); i++)
	{
		agents[i].act();
	}
	for (int i = 0; i < agents.size(); i++)
	{
		agents[i].update();
	}

	//	if (++nIteration == nIterTimes_)  break;
	//	// if(CompareState() > threshold ) break;
	//}

	std::cout << "evolve here ...\n";

	return 0;
}

int suStructrueOptimizer::run()
{
	int iteraTimes = 0;
	while (true) {
		evolve();
		/*outOofemfile("D:\\oofem\\build2.3\\Debug\\oofemOutFile.txt");
		useOofem();
		readFeedback();*/
		/*auto testIt = nodeArr_.begin();
		for (; testIt != nodeArr_.end(); testIt++) {
			if ((*testIt)->out == 1)
				std::cout << (*testIt)->morton << ":" << (*testIt)->strain << std::endl;
		}*/
	
		if (++iteraTimes == nIterTimes_) break;
	}
	export_stl_with_metaball("d:/test.stl", nodeArr_);





	return 0;
}

void suStructrueOptimizer::outOofemfile(std::string outAddress)
{
	auto outIt = nodeArr_.begin();
	int voxelNumber = 0;
	for (; outIt != nodeArr_.end(); outIt++) {
		if (((*outIt)->label_ == SU::BOUNDARY_CELL || (*outIt)->label_ == SU::BOUNDARY_CELL_SPECIAL 
			|| (*outIt)->label_ == SU::INTERIOR_CELL) && (*outIt)->out == true)
			voxelNumber++;
	}
	/*int count = 0;
	for (outIt = nodeArr_.begin(); outIt != nodeArr_.begin(); outIt++) {
		if ((*outIt)->label_ == SU::INTERIOR_CELL)
			count++;
	}
	std::cout << count << std::endl;
	system("pause");*/
	std::fstream outfile;
	outfile.open(outAddress, std::ios::out);
	outfile << "Majnun.out" << std::endl;
	outfile << "test of Brick elements with nlgeo 1(strain is the Green-Lagrangian strain) rotated as a rigid body" << std::endl;
	outfile << "#NonLinearStatic  nmsteps 1 nsteps 1 " << std::endl;
	outfile << "#LinearStatic  nmsteps 1 nsteps 1 " << std::endl;
	outfile << "#nsteps 5 rtolv 1.e-6 stiffMode 1 controlmode 1 maxiter 100" << std::endl;
	outfile << "#vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << std::endl;
	outfile << "#domain 3d" << std::endl;
	outfile << "#OutputManager tstep_all dofman_all element_all" << std::endl;
	outfile << "LinearStatic nsteps 3 nmodules 1" << std::endl;
	outfile << "vtkxml tstep_all domain_all primvars 1 1 vars 2 4 1 stype 1" << std::endl;
	outfile << "domain 3d" << std::endl;
	outfile << "OutputManager tstep_all dofman_all element_all" << std::endl;
	outfile << "ndofman " << pow((pow(2, nLevel_) + 1), 3) << " nelem " << voxelNumber << " ncrosssect 1 nmat 1 nbc 2 nic 0 nltf 1 " << std::endl;
	outfile.close();
	float x = (*pVolume_).bbMax_.data()[0] - (*pVolume_).bbMin_.data()[0];
	float y = (*pVolume_).bbMax_.data()[1] - (*pVolume_).bbMin_.data()[1];
	float z = (*pVolume_).bbMax_.data()[2] - (*pVolume_).bbMin_.data()[2];
	outCoor(x, y, z, nLevel_, outAddress);

	int number = 0;
	outIt = nodeArr_.begin();
	for (; outIt != nodeArr_.end(); outIt++)
	{
		if (((*outIt)->label_ == SU::BOUNDARY_CELL || (*outIt)->label_ == SU::BOUNDARY_CELL_SPECIAL 
			|| (*outIt)->label_ == SU::INTERIOR_CELL) && (*outIt)->out == true)
		{
			std::fstream outfile;
			outfile.open(outAddress, std::ios::app);
			int x, y, z;
			suMorton::decode(x, y, z, (*outIt)->morton, nLevel_);
			//outfile << pNode->xLocCode_ << "  " << pNode->yLocCode_ << "  " << pNode->zLocCode_ << "  " << pNode->level_ << "  " << pNode->label_ << "  ";
			//outfile << IndexX << "  " << IndexY << "  " << IndexZ;
			outfile << "LSpace " << ++number << "	 nodes  8 ";
			SU::voxel_output* asd = new SU::voxel_output(x, y, z, nLevel_,outAddress);
			outfile.close();
			//voxel_output asd(pChildNode, level);
			asd->output_point1();
			delete asd;
		}
	}

	outfile.open(outAddress, std::ios::app);
	outfile << "SimpleCS 1" << std::endl << "IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0" << std::endl << "BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.0"
		<< std::endl << "BoundaryCondition  2 loadTimeFunction 1 prescribedvalue 0.5" << std::endl << "PiecewiseLinFunction 1 npoints 2 t 2 0. 1000. f(t) 2 0. 1000." << std::endl;
	outfile.close();
}

void suStructrueOptimizer::useOofem()
{
	system("cd /d D:\\oofem\\build2.3\\Debug &oofem.exe&oofem -f oofemOutFile.txt");
}

void suStructrueOptimizer::outCoor(float maxX, float maxY, float maxZ, int level, std::string coorAddress)
{
	double lines = pow(2, level);//计算每一行有多少体素   dxdydz分别是三个方向上每个体素的尺寸
	double dx = maxX / lines;
	double dy = maxY / lines;
	double dz = maxZ / lines;

	std::fstream outfile;
	outfile.open(coorAddress, std::ios::app);
	int number = 1;
	lines += 1;
	for (int i = 0; i < lines; i++)//输出每个节点的坐标
	{
		for (int j = 0; j < lines; j++)
		{
			for (int k = 0; k < lines; k++)
			{
				if (k == 0)
				{
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << "       bc 3 2 2 2" << std::endl;
				}
				else
					outfile << "node   " << number++ << "   coords 3     " << j *dx << "  " << k *dy << "  " << i *dz << std::endl;
			}
		}
	}
	outfile.close();
}

void suStructrueOptimizer::readFeedback()
{
	read_point_information("D:\\oofem\\build2.3\\Debug\\Majnun.out.m0.1.vtu");
	float x = (*pVolume_).bbMax_.data()[0] - (*pVolume_).bbMin_.data()[0];
	float y = (*pVolume_).bbMax_.data()[1] - (*pVolume_).bbMin_.data()[1];
	float z = (*pVolume_).bbMax_.data()[2] - (*pVolume_).bbMin_.data()[2];
	read_point_coor("D:\\oofem\\build2.3\\Debug\\Majnun.out.m0.1.vtu",x,y,z,nLevel_);
	assignment(nodeArr_);
}

void suStructrueOptimizer::read_point_information(std::string address)
{
	all_point_mises_strain.clear();
	std::fstream in;
	char read_temp;
	in.open(address, std::ios::in);
	while (!in.eof())
	{
		in >> read_temp;
		if (read_temp == 'n')
		{
			in >> read_temp;
			if (read_temp == 'e')
			{
				in >> read_temp;
				if (read_temp == 'n')
				{
					in >> read_temp;
					if (read_temp == 't')
					{
						in >> read_temp;
						if (read_temp == 's')
						{
							in >> read_temp;
							if (read_temp == '=')
							{
								in >> read_temp;
								if (read_temp == '\"')
								{
									in >> read_temp;
									if (read_temp == '9')
									{
										for (;;)
										{
											in >> read_temp;
											//if (read_temp == '<')
											//break;
											int break_ = 0;
											if (read_temp == '>')
											{
												point_strain inf_temp;
												double temp_vector[9];
												char a[12];
												int count = 0;
												for (;;)
												{
													//in >> read_temp;
													in >> read_temp;
													if (read_temp == '<')
													{
														break_++;
														break;
													}
													if (read_temp == '-')
													{
														for (int i = 0; i < 12; i++)
														{
															in >> a[i];
														}
														temp_vector[count] = -1 * trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]);
													}
													else
													{
														a[0] = read_temp;
														for (int i = 1; i < 12; i++)
														{
															in >> a[i];
														}
														temp_vector[count] = trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]);
													}
													count++;
													if (count == 9)
													{
														count = 0;
														/*inf_temp.x_strain = pow(temp_vector[0] * temp_vector[0] + temp_vector[1] * temp_vector[1] + temp_vector[2] * temp_vector[2], 0.5);
														inf_temp.y_strain = pow(temp_vector[3] * temp_vector[3] + temp_vector[4] * temp_vector[4] + temp_vector[5] * temp_vector[5], 0.5);
														inf_temp.z_strain = pow(temp_vector[6] * temp_vector[6] + temp_vector[7] * temp_vector[7] + temp_vector[8] * temp_vector[8], 0.5);
														float strain_temp = pow(0.5*(pow(inf_temp.x_strain - inf_temp.y_strain, 2) + pow(inf_temp.y_strain - inf_temp.z_strain, 2) +
														pow(inf_temp.z_strain-inf_temp.x_strain, 2)), 0.5);*/
														float strain_temp = pow(pow(temp_vector[0] + temp_vector[3] + temp_vector[6], 2)
															+ pow(temp_vector[1] + temp_vector[4] + temp_vector[7], 2) + pow(temp_vector[2] + temp_vector[5] + temp_vector[8], 2), 0.5);
														std::cout << strain_temp << " ";
														all_point_mises_strain.push_back(strain_temp);
													}
												}
												break;
											}
											if (break_)
												break;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	in.close();
}

void suStructrueOptimizer::read_point_coor(std::string address, float box_maxx, float box_maxy, float box_maxz, int box_level)
{
	point_morton.clear();
	{
		float dx_ = box_maxx / pow(2, box_level);
		float dy_ = box_maxy / pow(2, box_level);
		float dz_ = box_maxz / pow(2, box_level);
		std::fstream in;
		char read_temp;
		in.open(address, std::ios::in);
		while (!in.eof())
		{
			in >> read_temp;
			//cout << read_temp;
			if (read_temp == 'i')
			{
				in >> read_temp;

				if (read_temp == 'i')
				{
					in >> read_temp;
					if (read_temp == '"')
					{
						in >> read_temp;
						//cout << read_temp;
						if (read_temp == '>')
						{
							float point_inf_temp[3];
							char a[12];
							int count = 0;
							for (;;)
							{
								//in >> read_temp;
								in >> read_temp;
								if (read_temp == '<')
								{
									break;
								}
								if (read_temp == '-')
								{
									for (int i = 0; i < 12; i++)
									{
										in >> a[i];
									}
									point_inf_temp[count] = -1 * trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]);
								}
								else
								{
									a[0] = read_temp;
									for (int i = 1; i < 12; i++)
									{
										in >> a[i];
									}
									point_inf_temp[count] = trans(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11]);
								}
								count++;
								if (count == 3)
								{
									count = 0;
									int tran[3];
									tran[0] = (point_inf_temp[0] + 0.5*dx_) / dx_;
									tran[1] = (point_inf_temp[1] + 0.5*dy_) / dy_;
									tran[2] = (point_inf_temp[2] + 0.5*dz_) / dz_;
									point_morton.push_back(suMorton::encode(tran[0], tran[1], tran[2], box_level + 1));
									std::cout << tran[0] << " " << tran[1] << " " << tran[2] << " " 
										<< suMorton::encode(tran[0], tran[1], tran[2], box_level + 1) << std::endl;

									//cout << inf_temp.x_strain << " " << inf_temp.y_strain << " " << inf_temp.z_strain << endl;
								}
							}
							break;
						}
					}
				}
			}
		}
		in.close();
		ass_point();
	}
}

double trans(char a1, char a2, char a3, char a4, char a5, char a6, char a7, char a8, char a9, char a10, char a11, char a12)
{
	double a = 0;
	a += asc2(a1);
	a += 0.1*asc2(a3);
	a += 0.01*asc2(a4);
	a += 0.001*asc2(a5);
	a += 0.0001*asc2(a6);
	a += 0.00001*asc2(a7);
	a += 0.000001*asc2(a8);
	int b = 0;
	b += asc2(a11) * 10 + asc2(a12);
	if (a10 == '-')
	{
		b = b*(-1);
	}
	a = a*pow(10, b);
	return a;
}

int asc2(char a)
{
	switch (a)
	{
	case '0':
		return 0;
		break;
	case '1':
		return 1;
		break;
	case '2':
		return 2;
		break;
	case '3':
		return 3;
		break;
	case '4':
		return 4;
		break;
	case '5':
		return 5;
		break;
	case '6':
		return 6;
		break;
	case '7':
		return 7;
		break;
	case '8':
		return 8;
		break;
	case '9':
		return 9;
		break;
	}
}

void suStructrueOptimizer::ass_point()
{
	for (int i = 0; i < point_morton.size(); i++)
	{
		point_mises_strain temp;
		temp.mises_strain = all_point_mises_strain[i];
		temp.point_morton = point_morton[i];
		sort_vector.push_back(temp);
	}
}

void suStructrueOptimizer::assignment(std::vector<SU::OctNode*>& assVector)
{
	auto assIt = assVector.begin();
	for (; assIt != assVector.end(); assIt++)//对所有的auto_cell赋应变值
	{
		(*assIt)->strain = return_max_strain((*assIt)->morton, nLevel_);
		//std::cout << auto_cell[ca[i]].strain << "  ";
		//std::cout << assIt->strain << ' ';
	}
}

float suStructrueOptimizer::return_max_strain(int morton_code, int level)
{
	int x, y, z;
	suMorton::decode(x,y,z,morton_code,level);
	//先计算对应的八个顶点的莫顿序
	int point_morton_[8];
	point_morton_[0] = suMorton::encode(x, y, z, level + 1);
	point_morton_[1] = suMorton::encode(x + 1, y, z, level + 1);
	point_morton_[2] = suMorton::encode(x, y + 1, z, level + 1);
	point_morton_[3] = suMorton::encode(x, y, z + 1, level + 1);
	point_morton_[4] = suMorton::encode(x + 1, y + 1, z + 1, level + 1);
	point_morton_[5] = suMorton::encode(x, y + 1, z + 1, level + 1);
	point_morton_[6] = suMorton::encode(x + 1, y, z + 1, level + 1);
	point_morton_[7] = suMorton::encode(x + 1, y + 1, z, level + 1);
	float max_strain = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < point_morton.size(); j++)
		{
			//cout << sort_vector[j].point_morton << " " << point_morton[i]<<":";
			if (sort_vector[j].point_morton == point_morton_[i])
			{
				max_strain = max_strain < abs(sort_vector[j].mises_strain) ? abs(sort_vector[j].mises_strain) : max_strain;
				//cout << "qwe"<<sort_vector[j].mises_strain << endl;
				break;
			}
		}
	}
	return max_strain;
}


bool suStructrueOptimizer::export_stl_with_metaball(const char * fileName, std::vector<SU::OctNode*>& stlVector)
{
	if (stlVector.size() < 2)
	{
		std::cout << "Open a model and make voxelization first!\n";
		return false;
	}
	//generate bounding box
	std::vector<Eigen::Vector3f> bbox;

	Eigen::Vector3d maxP, minP;
	maxP << pVolume_->bbMax_[0], 
		pVolume_->bbMax_[1],
		pVolume_->bbMax_[2];
	minP << pVolume_->bbMin_[0] ,
		pVolume_->bbMin_[1] ,
		pVolume_->bbMin_[2] ;

	Eigen::Vector3f fminP, fmaxP;
	fminP << minP(0), minP(1), minP(2);
	fmaxP << maxP(0), maxP(1), maxP(2);

	int nVolume_per_dim = pow(2, nLevel_);
	float volume_size = (maxP - minP).maxCoeff() / nVolume_per_dim;
	//volume_size = 0.0f;
	//extrend bbox with one voxel 
	Eigen::Vector3f voxSize;
	voxSize << volume_size, volume_size, volume_size;
	voxSize *= 2;

	bbox.push_back(fminP);// - voxSize);
	bbox.push_back(fmaxP);// +voxSize);

						  //generate metaballs on the model surface.
	std::vector<SU::METABALL> mballs;
	for (unsigned int i = 0; i <stlVector.size(); i++)
	{
		SU::METABALL m;

		//if (stlVector[i]->label_ == SU::BOUNDARY_CELL || stlVector[i]->label_ == SU::BOUNDARY_CELL_SPECIAL)
		if (/*!(stlVector[i]->out) &&*/ (stlVector[i]->out == false)  && stlVector[i]->label_ == SU::INTERIOR_CELL)//mellballs only include the voxels been cut
		{
			/*SU::Point cP = (stlVector[i]->min_ + stlVector[i]->max_) / 2;
			//TODO cP的获得方法要改掉
			int x, y, z;
			suMorton::decode(x, y, z, stlVector[i]->morton, nLevel_);
			cP.x = pVolume_->bbMin_.data()[0] + (pVolume_->bbMax_.data()[0] - pVolume_->bbMin_.data()[0]) / pow(2, nLevel_)*x;
			cP.y = pVolume_->bbMin_.data()[1] + (pVolume_->bbMax_.data()[1] - pVolume_->bbMin_.data()[1]) / pow(2, nLevel_)*y;
			cP.z = pVolume_->bbMin_.data()[2] + (pVolume_->bbMax_.data()[2] - pVolume_->bbMin_.data()[2]) / pow(2, nLevel_)*z;
			*/
			SU::Point cP = stlVector[i]->center();
			m.position << cP.x, cP.y, cP.z;
			m.squaredRadius = pow(volume_size / 2, 2);     //  metaball.r^2 = (volume box size / 2)^2
														   //std::cout << "r^2 = " << m.squaredRadius << std::endl;

			mballs.push_back(m);
		}


	}

	///resample
	/*int nSample = 5000;

	if (nSample > mballs.size()) nSample = mballs.size();

	std::vector<int> idxArr;
	std::vector<SU::METABALL> samples;
	for (int i = 0; i <mballs.size(); i++) idxArr.push_back(i);

	for (int i = 0; i < nSample; ++i) {
	int index = rand() % (mballs.size() - i) + i;
	std::swap(idxArr[i], idxArr[index]);
	samples.push_back(mballs[index] );
	}	*/

	/*
	SU::METABALL metaball;
	metaball.position = (bbox[0] + bbox[1]) / 2;
	metaball.squaredRadius = pow(volume_size / 2, 2);
	mballs.push_back(metaball);*/

	int resolution = nVolume_per_dim * 5;
	std::cout << "Num of Metaballs: " << mballs.size() << std::endl;
	std::cout << "Dim of Grid: " << resolution << "^ 3" << std::endl;
	std::cout << "stlVector.size()=" << stlVector.size() << std::endl;

	
	//取势场中threshold=1的等值面	
	if (mballs.size())
		SU::write_metaball_to_stl(fileName, /*samples*/mballs, fThresholdMC, resolution, bbox,pVolume_->mesh_);

	return true;
}