
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_GENERAL_SECTION_CXX

#include "ClassGeneralSection.hxx"

namespace CoefficientRepartition
{
  // Constructors.
  ClassGeneralSection::ClassGeneralSection()
    : ClassCoefficientRepartitionBase(), size_bin_(0), composition_bin_(0)
  {
    return;
  }


  ClassGeneralSection::ClassGeneralSection(const int &size_bin, const int &composition_bin)
    : ClassCoefficientRepartitionBase(), size_bin_(size_bin), composition_bin_(composition_bin)
  {
    if (size_bin_ < 0 || size_bin_ >= Nb_)
      throw CoefficientRepartition::Error("Out of range size bin.");

    if (composition_bin_ < 0 || composition_bin_ >= Nc_)
      throw CoefficientRepartition::Error("Out of range composition bin.");

    return;
  }


  // Destructor.
  ClassGeneralSection::~ClassGeneralSection()
  {
    return;
  }


  // Init.
  void ClassGeneralSection::Init(const string type)
  {
    Ops::Ops ops(configuration_file_);

    ops.SetPrefix(type + ".diameter.");

    //Number of fractions
    int Nfrac_;
    
    //Number of species
    int Nspecies_;
    Nspecies_ =ClassCoefficientRepartitionBase::GetNs();
    // Get diameter bounds.
    vector<double> diameter;
    if (ops.Exists("bound"))
      {
        ops.Set("bound", "", diameter);
      }
    else
      {
	Nfrac_ = ops.Get<int>("Nf");
        Nb_ = ops.Get<int>("Nb");
        double dmin = ops.Get<double>("min");
        double dmax = ops.Get<double>("max");
        diameter.resize(Nb_ + 1);
        double q = dmax / dmin;
        for (int j = 0; j < Nb_ + 1; j++)
          diameter[j] = dmin * pow(q, double(j) / double(Nb_));
      }

    if(Nb_==7)
    {
      diameter[0]=0.001;
      diameter[1]=0.00316;
      diameter[2]=0.01;
      diameter[3]=0.0398;
      diameter[4]=0.1585;
      diameter[5]=0.631;
      diameter[6]=2.5119;
      diameter[7]=10.0;
    }
    else
    {
      if(Nb_==5)
      {
      diameter[0]=0.01;
      diameter[1]=0.0398;
      diameter[2]=0.1585;
      diameter[3]=0.631;
      diameter[4]=2.5119;
      diameter[5]=10.0;
	  // diameter[0]=0.01;
	  // diameter[1]=0.1585;
	  // diameter[2]=0.4;
	  // diameter[3]=1.0;
	  // diameter[4]=2.5119;
	  // diameter[5]=10.0;
      }
      else
	if(Nb_==6)
	{
	  diameter[0]=0.01;
	  diameter[1]=0.0398;
	  diameter[2]=0.1585;
	  diameter[3]=0.4;
	  diameter[4]=1.0;
	  diameter[5]=2.5119;
	  diameter[6]=10.0;
	}
    }
    Nb_ = int(diameter.size()) - 1;//in case of diameter from existing files

    if (! Nb_ > 0)
      throw CoefficientRepartition::Error("Number of size sections should be strictly positive.");    

    vector<double> mass(Nb_ + 1);
    for (int j = 0; j < Nb_ + 1; j++)
      mass[j] = ClassCoefficientRepartitionBase::density_ * PI_6 * diameter[j] * diameter[j] * diameter[j];

    diameter_.Reallocate(Nb_);
    mass_.Reallocate(Nb_);
    for (int j = 0; j < Nb_; j++)
      for (int k = 0; k < 2; k++)
        {
          diameter_(j).PushBack(diameter[j + k]);
          mass_(j).PushBack(mass[j + k]);
        }

   // bool is_auto = ops.Get<bool>("section_definition_auto", "", true);

   if(Nfrac_!=0)//calculate composition sections automatically!!ops.Exists("autofraction")
     {
       //Get auto fraction number for every species
     
       cout<<"Auto composition sections:"<<endl;
       cout<<"Nfrac_= "<<Nfrac_<<"  Nspecies_= "<<Nspecies_<<endl;
       double sumfrac = 0;
       Nc_ = 0;
       int setb=0;
       vector<int> counter(Nspecies_ - 1);
       vector<double> frac_bound(Nfrac_ + 1);//fraction bounds for primary
       vector<double> compositions;
       //auto define the fraction bounds
       if (Nfrac_ == 3)
         {
           // frac_bound[0]=0.0;
           // frac_bound[1]=0.1;
           // frac_bound[2]=0.9;
           // frac_bound[3]=1.0;
           frac_bound[0]=0.0;
           frac_bound[1]=0.2;
           frac_bound[2]=0.8;
           frac_bound[3]=1.0;
         }
       else if(Nfrac_ == 4)
         {
           frac_bound[0]=0.0;
           frac_bound[1]=0.1;
           frac_bound[2]=0.5;
           frac_bound[3]=0.9;
           frac_bound[4]=1.0;
         }
       else
         {
           for(int i=0;i<Nfrac_ +1 ;i++)
             frac_bound[i] = ((double)i)/((double)Nfrac_);
         }

       //calculate the maximum fraction combinations
       for(int i=0; i<Nfrac_; i++)
         {
           for(int s=0; s<Nspecies_ - 1; s++)
             counter[s] = 0;//initial the counter
           //when the index counter of second species reaches its top, 
           //move to the next fraction bin of first species
           if(Nspecies_>2)
             {
               while(counter[1] <= Nfrac_ -1)
                 {//take the base fraction bounds of current bin of first species
                   sumfrac=frac_bound[i];
                   for(int s=1; s<Nspecies_-1;s++)
                     {
                       int j=counter[s];//the fraction bin index for species s
                       sumfrac=sumfrac+frac_bound[j];//calculate one possible combination
                     }
                   if(sumfrac<1.0)
                     {
                       Nc_++;
                       //for first group
                       compositions.push_back(frac_bound[i]);
                       compositions.push_back(frac_bound[i+1]);
                       for(int g=1; g<Nspecies_-1;g++)
                         {//save possible combinations
                           int j=counter[g];//the fraction list index for group g
                           compositions.push_back(frac_bound[j]);
                           compositions.push_back(frac_bound[j+1]);
                         }
                       //last group with default fraction section [0,1]
                       compositions.push_back(0.0);
                       compositions.push_back(1.0);	      
                     }
                   //when the second last species hasn't reaches its top, 
                   if(counter[Nspecies_-2]<=Nfrac_)
                     counter[Nspecies_-2]++;//move the index of second last species
                   for(int g=3;g<Nspecies_;g++)
                     {//check every neighbor counter,[2,Ngroup_aer-2] form back to forward
                       int j=Nspecies_+1-g;
                       sumfrac=frac_bound[counter[j-1]]+frac_bound[counter[j]];
                       //the bottom sum of two neighbor group fraction is already too big
                       if(sumfrac>=1.0)
                         {
                           //reset all j following counter
                           for(int s=j;s<Nspecies_-1;s++)
                             counter[s]=0;
                           //increase the
                           counter[j-1]++;
                         }
                     }
                 }
             }
           else
             {//n case of only two
               if(Nspecies_==2)
                 {
                   Nc_=Nfrac_;
                   //for first group
                   compositions.push_back(frac_bound[i]);
                   compositions.push_back(frac_bound[i+1]);
                   //for second group
                   compositions.push_back(0.0);
                   compositions.push_back(1.0);
                 }
               else
                 {//one species
                   Nc_=1;
                   Nfrac_=1;
                   compositions.push_back(0.0);
                   compositions.push_back(1.0);	    
                 }
             }
         }
       cout<<"Number of composition sections= "<<Nc_<<endl;
       fraction_index_.Reallocate(Nc_);
       fraction_.Reallocate(Nc_);
       int iter=0;
       for(int c=0; c<Nc_; c++)
         {
           fraction_index_(c).Reallocate(Nspecies_);
           fraction_(c).Reallocate(Nspecies_);
           for(int s=0;s<Nspecies_;s++)
             for(int k=0; k<2; k++)
               {
                 cout<<"fraction_("<<c<<","<<s<<")="<<compositions[iter]<<endl;
                 fraction_(c)(s).PushBack(compositions[iter]);//c is the index of composition_section
                 fraction_index_(c)(s).PushBack(s);
                 iter++;
               }
         }
     }
   else
     {

      // Get composition sections.
      ops.SetPrefix(type + ".");

      // Number of composition sections.
      Nc_ = 0;
      while (ops.Exists("composition_section_" + to_str(Nc_)))
	Nc_++;

      if (! Nc_ > 0)
	throw CoefficientRepartition::Error("Number of composition sections should be strictly positive.");    

      fraction_index_.Reallocate(Nc_);
      fraction_.Reallocate(Nc_);

      // Read definition of each composition section.
      for (int c = 0; c < Nc_; c++)
	{
	  ops.SetPrefix(type + ".");

	  vector<string> vtmp = ops.GetEntryList("composition_section_" + to_str(c));

	  fraction_index_(c).Reallocate(int(vtmp.size()));//vtmp.size()=number groups
	  fraction_(c).Reallocate(int(vtmp.size()));

	  ops.SetPrefix(type + ".composition_section_" + to_str(c) + ".");
	  for (int i = 0; i < int(vtmp.size()); i++)
	    {
	      vector<double> fraction;
	      ops.Set(vtmp[i], "", fraction);
	      fraction_(c)(i).PushBack(fraction[0]);//c is the index of composition_section
	      fraction_(c)(i).PushBack(fraction[1]);//i is the index of species
	      vector<string> ls = split(vtmp[i].substr(1), "+");
	      fraction_index_(c)(i).Reallocate(int(ls.size()));//ls.size()=number of species within each groups
	      for (int j = 0; j < int(ls.size()); j++)//j is always unkown!!
	      {
		fraction_index_(c)(i)(j) = convert<int>(ls[j]);
	      }
	    }
	}
    }
    ops.Close();
  }


  // Copy.
  void ClassGeneralSection::Copy(const ClassGeneralSection &gs)
  {
    size_bin_ = gs.size_bin_;
    composition_bin_ = gs.composition_bin_;
  }
  
  // Get methods.
  int ClassGeneralSection::GetNd() const
  {
    return fraction_(composition_bin_).GetSize();
  }


  int ClassGeneralSection::GetSizeBin() const
  {
    return size_bin_;
  }


  int ClassGeneralSection::GetCompositionBin() const
  {
    return composition_bin_;
  }


  double ClassGeneralSection::GetMassBoundary(const int &i) const
  {
    return mass_(size_bin_)(i);
  }


  double ClassGeneralSection::GetDiameterBoundary(const int &i) const
  {
    return diameter_(size_bin_)(i);
  }


  void ClassGeneralSection::CollectFractionIndexBoundary(const int &i, Vector1I &fraction_index) const
  {
    fraction_index = fraction_index_(composition_bin_)(i);
  }


  void ClassGeneralSection::CollectFractionBoundary(const int &i, Vector1T &fraction) const
  {
    fraction = fraction_(size_bin_)(i);
  }


  int ClassGeneralSection::GetNdStatic(const int &c)
  {
    return fraction_(c).GetSize();
  }


  int ClassGeneralSection::GetNb()
  {
    return Nb_;
  }


  int ClassGeneralSection::GetNc()
  {
    return Nc_;
  }


  void ClassGeneralSection::CollectDiameter(Vector1T &diameter)
  {
    diameter.Reallocate(Nb_ + 1);
    for (int b = 0; b < Nb_; b++)
      diameter(b) = diameter_(b)(0);
    diameter(Nb_) = diameter_(Nb_ - 1)(1);
  }


  void ClassGeneralSection::CollectMass(Vector1T &mass)
  {
    mass.Reallocate(Nb_ + 1);
    for (int b = 0; b < Nb_; b++)
      mass(b) = mass_(b)(0);
    mass(Nb_) = mass_(Nb_ - 1)(1);
  }


  void ClassGeneralSection::CollectFractionIndexBoundaryStatic(const int &c, const int &i, Vector1I &fraction_index)
  {
    fraction_index = fraction_index_(c)(i);
  }


  void ClassGeneralSection::CollectFractionBoundaryStatic(const int &c, const int &i, Vector1T &fraction)
  {
    fraction = fraction_(c)(i);
  }


  // Randomly generate particles within section.
  void ClassGeneralSection::GenerateRandomParticle(ClassParticle &p) const
  {
    return generate_random_particle(p);
  }


  void ClassGeneralSection::generate_random_particle(ClassParticle &p) const
  {//fund the worst problem of negtive fraction particles
    Vector1T &mass = mass_(size_bin_);//size_bin_ define current size bin index
    double frac_limt(double(0));
    double frac_avalibe(double(1));
    double frac_sum(double(0));
    p.mass_ = mass(0) + drand48() * (mass(1) - mass(0));//particle total mass
    Vector<bool> has_fraction(this->Ns_);//array to verify correspondence of each species with their fraction section
    has_fraction.Fill(false);
    Vector2I &fraction_index = fraction_index_(composition_bin_);//composition_bin_ define current composition ID (s)(0)=s
    Vector2T &fraction = fraction_(composition_bin_);//fraction bounds of each compositions (s)(0)/(1)
    p.fraction_.Zero();
    int count_has_fraction(0);
    //fraction_index.GetSize()=Nspecies_-1?N_groups-1
    for (int i = 0; i < fraction_index.GetSize(); i++)
      {
	frac_limt+=fraction(i)(0);
      }
      frac_avalibe=(double(1))-frac_limt;
      frac_limt=frac_avalibe;
    //Generate fraction before last species
    for (int i = 0; i < fraction_index.GetSize()-1; i++)
      {
	double frac;
	double frac_section=fraction(i)(1) - fraction(i)(0);
	double real_avalible=(frac_section<frac_avalibe)?frac_section:frac_avalibe;//mini value

	frac= fraction(i)(0) + drand48() * real_avalible;
	frac_avalibe-=(frac-fraction(i)(0));
	
        double frac_rest(double(1));

        int s = fraction_index(i)(0);
        p.fraction_(s) = frac * frac_rest;
        frac_sum += frac;
      }
      
    if(frac_sum>1.0)
    {
      cout<<"fraction error1: frac_sum="<<frac_sum<<endl;
      for(int i = 0; i < fraction_index.GetSize(); i++)
      {
	int s = fraction_index(i)(0);
	cout<<"frac_"<<s<<"="<<p.fraction_(s)<<endl;
	cout<<"("<<fraction(s)(0)<<","<<fraction(s)(1)<<")"<<endl;
      }
     exit(2);
    }

    //Generate fraction for last species Ns_
   int s=fraction_index.GetSize()-1;
   p.fraction_(s)=1.0-frac_sum;
   frac_sum+=p.fraction_(s);

   //check mass and fraction for generated particles
  for(int s=0;s<this->Ns_;s++)
  {
    if(p.fraction_(s)>=fraction(s)(0)&&p.fraction_(s)<=fraction(s)(1))
    {
      has_fraction(s) = true;
      count_has_fraction++;
    }
  }

     if (count_has_fraction < this->Ns_|| frac_sum!=1.0)
       {
	cout<<"fraction error2: frac_sum="<<frac_sum<<endl;
	for(int i = 0; i < fraction_index.GetSize(); i++)
	{
	  int s = fraction_index(i)(0);
	  cout<<"frac_s="<<s<<"="<<p.fraction_(s)<<endl;
	  cout<<"("<<fraction(s)(0)<<","<<fraction(s)(1)<<")"<<endl;
	}
	exit(2);
      }
      
  }


  // Whether one particle is inside general section.
  bool ClassGeneralSection::HasParticle(const ClassParticle &p) const
  {
    return has_particle(p);
  }


  inline bool ClassGeneralSection::has_particle(const ClassParticle &p) const
  {
    Vector1T &mass = mass_(size_bin_);

    if (size_bin_ == (Nb_-1))
      {
        // If on last bin, consider particle to be in section even if
        // its mass exceeds the upper bound, thus not checked hereafter.
        if (p.mass_ < mass(0))
          return false;
      }
    else
    {
      if (p.mass_ < mass(0) || p.mass_ >= mass(1))
        return false;
    }

    Vector2I &fraction_index = fraction_index_(composition_bin_);
    Vector2T &fraction = fraction_(composition_bin_);

    for (int i = 0; i < fraction_index.GetSize(); i++)
      {
        double frac_sum(0);
        for (int j = 0; j < 1; j++)//fraction_index(i).GetSize()
	{// j=0 always 1 species for each group
          frac_sum += p.fraction_(fraction_index(i)(j));
	}

	if(frac_sum<1.0)
	{
	  if (frac_sum < fraction(i)(0) || frac_sum >= fraction(i)(1))
	    return false;
	}
	else
	{
	  if(fraction(i)(1)!=1.0)
	    return false;
	}
      }
    return true;
  }
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_GENERAL_SECTION_CXX
#endif
