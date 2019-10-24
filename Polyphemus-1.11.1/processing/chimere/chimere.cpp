#include <iostream>
using namespace std;

#include "Verdandi.hxx"
#include "Chimere.cxx"

#include "BaseDriver.cxx"
#include "BaseOutputSaver.cxx"

using namespace Polyphemus;

int main(int argc, char** argv)
{
  TRY;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  BaseDriver<float, Chimere, BaseOutputSaver<float, Chimere> >
    driver(argv[1]);

  driver.Run();

  END;

  return 0;
}
