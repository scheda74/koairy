
#ifndef COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_CXX

#include "ClassCoefficientRepartition.hxx"

namespace CoefficientRepartition
{
  // Constructors.
  ClassCoefficientRepartition::ClassCoefficientRepartition(const string type,
                                                           const int Nmc)
    : ClassCoefficientRepartitionBase(), Ncompute_(0), Nmc_(Nmc), Nsize_(0)
  {
    Ops::Ops ops(ClassCoefficientRepartitionBase::configuration_file_);
    ops.SetPrefix(type + ".");

    Nmc_ = ops.Get<int>("Number_monte_carlo", "", Nmc);

    if (! Nmc_ > 0)
      throw CoefficientRepartition::Error("Number of Monte Carlo should be strictly positive.");

    bool is_auto = ops.Get<bool>("section_definition_auto", "", true);

    if (is_auto)
      Nsize_ = ClassGeneralSection::GetNc() * ClassGeneralSection::GetNb();
    else
      while (ops.Exists("general_section_" + to_str(Nsize_)))
        Nsize_++;

    if (! Nsize_ > 0)
      throw CoefficientRepartition::Error("Number of general sections should be strictly positive.");

    general_section_.Reallocate(Nsize_);

    if (is_auto)
      for (int i = 0; i < Nsize_; i++)
        {

          int size_bin = i / ClassGeneralSection::GetNc();
          int composition_bin = i % ClassGeneralSection::GetNc();
          //cout<<i<<" size_bin="<<size_bin<<" composition_bin="<<composition_bin<<endl;
          ClassGeneralSection gs(size_bin, composition_bin);
          general_section_(i).Copy(gs);
        }
    else
      for (int i = 0; i < Nsize_; i++)
        {
          vector<int> itmp;
          ops.Set("general_section_" + to_str(i), "", itmp);
          int size_bin = itmp[0];
          int composition_bin = itmp[1];

          ClassGeneralSection gs(size_bin, composition_bin);
          general_section_(i).Copy(gs);
        }

    ops.Close();

    index_first_.Reallocate(Nsize_);
    index_second_.Reallocate(Nsize_);
    coefficient_.Reallocate(Nsize_);

    return;
  }


  // Destructor.
  ClassCoefficientRepartition::~ClassCoefficientRepartition()
  {
    return;
  }


  // Get methods.
  int ClassCoefficientRepartition::GetNcompute() const
  {
    return Ncompute_;
  }


  int ClassCoefficientRepartition::GetNsize() const
  {
    return Nsize_;
  }


  int ClassCoefficientRepartition::GetNmc() const
  {
    return Nmc_;
  }


  ClassGeneralSection* ClassCoefficientRepartition::GetGeneralSection(const int &i)
  {
    return &general_section_(i);
  }


  void ClassCoefficientRepartition::CollectIndexFirst(const int &i, Vector1I &index) const
  {
    index = index_first_(i);
  }
  

  void ClassCoefficientRepartition::CollectIndexSecond(const int &i, Vector1I &index) const
  {
    index = index_second_(i);
  }


  void ClassCoefficientRepartition::CollectCoefficient(const int &i, Vector1T &coefficient) const
  {
    coefficient = coefficient_(i);
  }


  // Clear.
  void ClassCoefficientRepartition::Clear()
  {
    Ncompute_ = 0;
    index_first_.Clear();
    index_second_.Clear();
    coefficient_.Clear();
  }


  // Compute repartition coefficients between general sections i1 and i2.
  void ClassCoefficientRepartition::Compute(const int &i1, const int &i2)
  {
    Vector1I count(Nsize_);
    count.Zero();
    cout<<"i1= "<<i1<<" i2= "<<i2<<endl;
    // Monte Carlo loop.
    for (int u = 0; u < Nmc_; u++)
      {
        ClassParticle particle1, particle2;
        general_section_(i1).generate_random_particle(particle1);
        general_section_(i2).generate_random_particle(particle2);
	
        ClassParticle particle12 = particle1.coagulate(particle2);
	int j;
        for (int i = 0; i < Nsize_; i++)
	{
          if (general_section_(i).has_particle(particle12))
            {
              count(i)++;
              break;
            }
            j=i;
	}
	//what if i=Nsize and has_particle never be ture?
	if(j==Nsize_-1)
	{
	  cout<<" p12.m="<<particle12.GetMass()<<endl;
	  cout<<general_section_(i1).GetMassBoundary(0)<<" "<<general_section_(i1).GetMassBoundary(1)<<endl;
	  Vector1T fraction;
	  particle12.CollectFraction(fraction);
	  for (int g = 0; g <5; g++)
	    cout<<fraction(g)<<" "<<g<<endl;
	  cout<<"error, no match "<<i1<<" "<<i2<<endl;
	  exit(2);
	}
      }

    for (int i = 0; i < Nsize_; i++)
    {
      if (count(i) > 0)
        {
          index_first_(i).PushBack(i1);
          index_second_(i).PushBack(i2);
          coefficient_(i).PushBack(double(count(i)) / double(Nmc_));
        }
        int Ncoef = coefficient_(i).GetSize();
	if(Ncoef==0)
	{
	  index_first_(i).PushBack(i);
          index_second_(i).PushBack(i);
          coefficient_(i).PushBack(double(0.0));
	 }
    }

    Ncompute_++;
  }


  void ClassCoefficientRepartition::ComputeAll()
  {
    int Ncouple_global = Nsize_ * (Nsize_ + 1) / 2;
    int Ncouple_local = Ncouple_global / this->Nrank_;
    int Ncouple_remaining = Ncouple_global % this->Nrank_;
    cout<<"Ncouple_global: "<<Ncouple_global<<endl;
    cout<<"Ncouple_local: "<<Ncouple_local<<endl;
    cout<<"Ncouple_remaining: "<<Ncouple_remaining<<endl;
    if (this->rank_ < Ncouple_remaining)
      Ncouple_local++;
    cout<<"Nrank_: "<<this->Nrank_<<endl;
    cout<<"Ncouple_local: "<<Ncouple_local<<endl;

    int Ncouple_offset = Ncouple_global / this->Nrank_ * this->rank_;
    cout<<"Ncouple_offset: "<<Ncouple_offset<<endl;

    for (int irank = 0; irank < this->rank_; irank++)
      if (irank < Ncouple_remaining)
        Ncouple_offset++;

    int i(0);
    Vector2I couple(Ncouple_global);
    for (int i1 = 0; i1 < Nsize_; i1++)
      for (int i2 = 0; i2 < i1 + 1; i2++)
        {
          couple(i).Reallocate(2);
          couple(i)(0) = i1;
          couple(i)(1) = i2;
          i++;
        }
    cout<<"i= "<<i<<endl;
    for (int i = 0; i < Ncouple_local; i++)
      {
        int j = Ncouple_offset + i;
        if (j != i)
          couple(i) = couple(j);
      }

    couple.Resize(Ncouple_local);

    cout << "Rank " << this->rank_ << "/" << this->Nrank_
         << " computes " << Ncouple_local << " couples from ("
         << couple(0)(0) << ", " << couple(0)(1) << ") to ("
         << couple(Ncouple_local - 1)(0) << ", " << couple(Ncouple_local - 1)(1) << ")" << endl;
   
    // Compute starting from last computed couple.
    for (int i = Ncompute_; i < Ncouple_local; i++)
    {
      cout<<" couples: "<<i<<" out of: "<<Ncouple_local<<"("<<couple(i)(0)<<","<<couple(i)(1)<<")"<<endl;
      Compute(couple(i)(0), couple(i)(1));
    }
   
  }


#ifdef WITH_MPI
  void ClassCoefficientRepartition::MPI_ComputeAll(const int recv_rank)
  {
    ComputeAll();

    MPI::COMM_WORLD.Barrier();

    for (int i = 0; i < Nsize_; i++)
      {
        Vector1I rcounts(this->Nrank_);
        rcounts.Zero();

        // Send size to receiver rank.
        rcounts(this->rank_) = coefficient_(i).GetSize();
        MPI::COMM_WORLD.Gather(rcounts.GetData() + this->rank_, 1, MPI::INT,
                               rcounts.GetData(), 1, MPI::INT, recv_rank);

        int rcount(0);
        Vector1I displs(this->Nrank_);
        displs.Zero();

        if (this->rank_ == recv_rank)
          {
            rcount = rcounts(0);
            for (int irank = 1; irank < this->Nrank_; irank++)
              {
                displs(irank) = rcount;
                rcount += rcounts(irank);
              }
          }

        MPI::COMM_WORLD.Bcast(&rcount, 1, MPI::INT, recv_rank);
        if (rcount == 0)
          continue;

        Vector1I recv_buf_int;
        Vector1T recv_buf_real;
        if (this->rank_ == recv_rank)
          {
            recv_buf_int.Reallocate(rcount);
            recv_buf_real.Reallocate(rcount);
          }

        int scount = coefficient_(i).GetSize();

        MPI::COMM_WORLD.Gatherv(index_first_(i).GetData(), scount, MPI::INT,
                                recv_buf_int.GetData(), rcounts.GetData(),
                                displs.GetData(), MPI::INT, recv_rank);
        if (this->rank_ == recv_rank)
          index_first_(i) = recv_buf_int;

        MPI::COMM_WORLD.Gatherv(index_second_(i).GetData(), scount, MPI::INT,
                                recv_buf_int.GetData(), rcounts.GetData(),
                                displs.GetData(), MPI::INT, recv_rank);
        if (this->rank_ == recv_rank)
          index_second_(i) = recv_buf_int;

        MPI::COMM_WORLD.Gatherv(coefficient_(i).GetData(), scount, COEFFICIENT_REPARTITION_MPI_REAL,
                                recv_buf_real.GetData(), rcounts.GetData(),
                                displs.GetData(), COEFFICIENT_REPARTITION_MPI_REAL, recv_rank);
        if (this->rank_ == recv_rank)
          coefficient_(i) = recv_buf_real;
      }

    MPI::COMM_WORLD.Barrier();
  }
#endif


#ifdef WITH_NETCDF
  // Read coefficients.
  void ClassCoefficientRepartition::ReadNetCDF(const string &input_file)
  {
    NcFile fnc(input_file.c_str(), NcFile::ReadOnly);

    if (! fnc.is_valid())
      throw CoefficientRepartition::Error("Invalid NetCDF file pointer.");

    // Clear data.
    Clear();

    // Check some dimensions.
    if (Nmc_ != int(fnc.get_dim("Nmc")->size()))
      throw CoefficientRepartition::Error("Values of Nmc differ between NetCDF file and current object.");

    if (this->Ns_ != int(fnc.get_dim("Ns")->size()))
      throw CoefficientRepartition::Error("Values of Ns differ between NetCDF file and current object.");

    if (Nsize_ != int(fnc.get_dim("Nsize")->size()))
      throw CoefficientRepartition::Error("Values of Nsize differ between NetCDF file and current object.");

    if (ClassGeneralSection::GetNb() != int(fnc.get_dim("Nb")->size()))
      throw CoefficientRepartition::Error("Values of Nb differ between NetCDF file and current object.");

    if (ClassGeneralSection::GetNc() != int(fnc.get_dim("Nc")->size()))
      throw CoefficientRepartition::Error("Values of Nc differ between NetCDF file and current object.");

    // Current number of couples computed.
    Ncompute_ = int(fnc.get_dim("Ncompute")->size());

    index_first_.Reallocate(Nsize_);
    index_second_.Reallocate(Nsize_);
    coefficient_.Reallocate(Nsize_);

    for (int i = 0; i < Nsize_; i++)
      {
        string dim_name("Ncoef_" + to_str(i));
        int Ncoef = int(fnc.get_dim(dim_name.c_str())->size());

        index_first_(i).Reallocate(Ncoef);
        index_second_(i).Reallocate(Ncoef);
        coefficient_(i).Reallocate(Ncoef);

        NcVar *var;
        string var_name;

        var_name = "index1_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(index_first_(i).GetData(), Ncoef);

        var_name = "index2_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(index_second_(i).GetData(), Ncoef);

        var_name = "coef_" + to_str(i);
        var = fnc.get_var(var_name.c_str());
        var->get(coefficient_(i).GetData(), Ncoef);
      }

    fnc.close();
  }

  void ClassCoefficientRepartition::WriteBIN(const string &output_file) const
  {  
    std::ofstream outputFile (output_file.c_str(), ios::out |  ios::binary);
    for (int i = 0; i < Nsize_; i++)
      {
	int Ncoef = coefficient_(i).GetSize();
	int ft_tag =4;
	outputFile.write((char*)&ft_tag, sizeof(int));
	outputFile.write((char*)&Ncoef, sizeof(int));
	outputFile.write((char*)&ft_tag, sizeof(int));
	
	for(int j=0;j<Ncoef;j++)
	{
	  int index1=index_first_(i).GetData()[j];
	  outputFile.write((char*)&ft_tag, sizeof(int));
	  outputFile.write((char*)&index1, sizeof(int));
	  outputFile.write((char*)&ft_tag, sizeof(int));
	  int index2=index_second_(i).GetData()[j];
	  outputFile.write((char*)&ft_tag, sizeof(int));
	  outputFile.write((char*)&index2, sizeof(int));
	  outputFile.write((char*)&ft_tag, sizeof(int));
	  double coeff=coefficient_(i).GetData()[j];
	  int coef_i=coeff*10000000;
	  outputFile.write((char*)&ft_tag, sizeof(int));
	  outputFile.write((char*)&coef_i, sizeof(int));
	  outputFile.write((char*)&ft_tag, sizeof(int));
	}
      }
      outputFile.close();
  }

  void ClassCoefficientRepartition::WriteTXT(const string &output_file) const
  {  
    std::ofstream outputFile (output_file.c_str(), ios::out );
    for (int i = 0; i < Nsize_; i++)
      {
	int Ncoef = coefficient_(i).GetSize();
	outputFile<<Ncoef<< endl;
	
	for(int j=0;j<Ncoef;j++)
	{
	  int index1=index_first_(i).GetData()[j];
	  outputFile<<index1<< endl;
	  int index2=index_second_(i).GetData()[j];
	  outputFile<<index2<< endl;
	  double coeff=coefficient_(i).GetData()[j];
	  outputFile<<coeff<< endl;
	}
      }
      outputFile.close();
  }

  // Write coefficients.
  void ClassCoefficientRepartition::WriteNetCDF(const string &output_file) const
  {
    NcFile fnc(output_file.c_str(), NcFile::New);

    NcDim *dim;

    dim = fnc.add_dim("Nmc", Nmc_);
    dim = fnc.add_dim("Ns", this->Ns_);
    dim = fnc.add_dim("Nsize", Nsize_);
    dim = fnc.add_dim("Nb", ClassGeneralSection::GetNb());
    dim = fnc.add_dim("Nc", ClassGeneralSection::GetNc());
    dim = fnc.add_dim("Ncompute", Ncompute_);

    for (int i = 0; i < Nsize_; i++)
      {
        int Ncoef = coefficient_(i).GetSize();//get the array dimension (three types of array have the same dimension)
        string dim_name("Ncoef_" + to_str(i));//the name of the dim is always changing
        dim = fnc.add_dim(dim_name.c_str(), Ncoef);//one dimension data

        NcVar *var;
        string var_name;

        var_name = "index1_" + to_str(i);
        var = fnc.add_var(var_name.c_str(), ncInt, dim);
        var->put(index_first_(i).GetData(), Ncoef);//writing data _ one dimension array

        var_name = "index2_" + to_str(i);
        var = fnc.add_var(var_name.c_str(), ncInt, dim);
        var->put(index_second_(i).GetData(), Ncoef);

        var_name = "coef_" + to_str(i);
        var = fnc.add_var(var_name.c_str(), COEFFICIENT_REPARTITION_NETCDF_REAL, dim);
        var->put(coefficient_(i).GetData(), Ncoef);
      }

    fnc.close();
  }
#endif
}

#define COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_CXX
#endif
