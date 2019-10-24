int i;

// Overrides if assigned or detected value is too large.
if (child_process_fork > MAX_CHILD_PROCESS)
  child_process_fork = MAX_CHILD_PROCESS;

int iquotient, iremainder;
iquotient = Model.GetNx() / child_process_fork;
iremainder = Model.GetNx() - child_process_fork * iquotient;

int iwidth[MAX_CHILD_PROCESS - 1];

for (i = 0; i < child_process_fork - 1; i++)
  {
    iwidth[i] = iquotient;
    if (iremainder > 0)
      {
        iwidth[i] += 1;
        iremainder -= 1;
      }
  }

ibeg[0] = 1;
for (i = 1; i < child_process_fork; i++)
  ibeg[i] = ibeg[i - 1] + iwidth[i - 1];

for (i = 0; i < child_process_fork - 1; i++)
  iend[i] = ibeg[i + 1] - 1;
iend[child_process_fork - 1] = Model.GetNx();

for (i = 0; i < child_process_fork; i++)
  {
    jbeg[i] = 1;
    jend[i] = Model.GetNy();
  }

// Memory for concentration fields (gaseous and aerosol), pH and
// in-cloud wet deposition fluxes (gaseous and aerosol).
for (i = 0; i < child_process_fork; i++)
  shared_memory_size[i] = (jend[i] - jbeg[i] + 1)
    * (iend[i] - ibeg[i] + 1)
    * (Model.GetNz()
       * (Nbin_aer * Ns_aer + Ns + 1)
       + (Nbin_aer * Ns_aer + Ns));

shared_memory.Init(child_process_fork, shared_memory_size);
