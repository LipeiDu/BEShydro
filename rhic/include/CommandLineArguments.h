//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef COMMANDLINEARGUMENTS_H_
#define COMMANDLINEARGUMENTS_H_

#include <stdbool.h>
#include <argp.h>

struct CommandLineArguments
{
  char *args[2]; /* ARG1 and ARG2 */
  bool runTest;
  bool runHydro;
  char *configDirectory; /* The -v flag */
  char *outputDirectory; /* Argument for -o */
};

error_t loadCommandLineArguments(int argc, char **argv, void * cli_params, const char *version, const char *address);

#endif /* COMMANDLINEARGUMENTS_H_ */
