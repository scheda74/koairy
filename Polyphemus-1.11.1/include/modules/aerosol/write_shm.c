int i1 = ibeg[c];
int i2 = iend[c];
int j1 = jbeg[c];
int j2 = jend[c];

int b, s, k, j, i, l;

double *data;

data = (double*) shmat(shared_memory.GetShmid(c), NULL, SHM_R | SHM_W);

l = 0;
for (i = i1 - 1; i < i2; i++)
  for (j = j1 - 1; j < j2; j++)
    for (k = 0; k < Nz; k++)
      for (s = 0; s < Ns_aer; s++)
        for (b = 0; b < Nbin_aer; b++)
          data[l++] = Concentration_aer(s, b, k, j, i);

for (i = i1 - 1; i < i2; i++)
  for (j = j1 - 1; j < j2; j++)
    for (k = 0; k < Nz; k++)
      for (s = 0; s < Ns; s++)
        data[l++] = Concentration(s, k, j, i);

shmdt(data);
