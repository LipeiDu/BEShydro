/*
 * CommandLineArguments.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "../include/CommandLineArguments.h"

const char *argp_program_version;
const char *argp_program_bug_address;

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
*/
static struct argp_option options[] =
{
		{"test",  't', "RUN_TEST", OPTION_ARG_OPTIONAL, "Run software tests"},
		{"hydro",  'h', "RUN_HYDRO", OPTION_ARG_OPTIONAL, "Run hydrodynamic simulation"},
		{"output",  'o', "OUTPUT_DIRECTORY", 0, "Path to output directory"},
		{"config", 'c', "CONFIG_DIRECTORY", 0, "Path to configuration directory"},
		{0}
};


/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt(int key, char *arg, struct argp_state *state) {
	struct CommandLineArguments *cli = (CommandLineArguments *)state->input;

	switch (key) {
	case 't':
		cli->runTest = true;
		break;
	case 'h':
		cli->runHydro = true;
		break;
	case 'o':
		cli->outputDirectory = arg;
		break;
	case 'c':
		cli->configDirectory = arg;
		break;
//	case ARGP_KEY_ARG:
//		if (state->arg_num >= 2) {
//			argp_usage(state);
//		}
//		arguments->args[state->arg_num] = arg;
//		break;
//	case ARGP_KEY_END:
//		if (state->arg_num < 2) {
//			argp_usage(state);
//		}
//		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
     that we accept.
*/
static char args_doc[] = "ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"Run -- A program to run a single viscous hydrodynamic simulation of a relativistic heavy ion collision";

/*
   The ARGP structure itself.
*/
static struct argp argp = {options, parse_opt, args_doc, doc};

/*
   The main function.
   Notice how now the only function call needed to process
   all command-line options and arguments nicely
   is argp_parse.
*/
error_t loadCommandLineArguments(int argc, char **argv, void * cli_params, const char *version, const char *address)
{
	struct CommandLineArguments * cli = (struct CommandLineArguments *)cli_params;

	argp_program_version = version;
	argp_program_bug_address = address;

  /* Set argument defaults */
	cli->runTest = false;
	cli->runHydro = false;
	cli->outputDirectory = NULL;
	cli->configDirectory = NULL;

  /* Where the magic happens */
  argp_parse (&argp, argc, argv, 0, 0, cli);

  return 0;
}
