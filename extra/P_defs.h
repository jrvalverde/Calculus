/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*   P_DEFS.H							*
*								*
*   Constant definitions for the Protein Identification and	*
*   Library Access Resource.					*
*								*
*								*
*   DESIGNED BY:						*
*	José Ramón Valverde Carrillo.				*
*								*
*   HISTORY:							*
*	J. R. Valverde (1989):					*
*	    Updated as needed while developing the project.	*
*	J. R. Valverde ( 15 - aug - 1989):			*
*	    Version 1 concluded.				*
*	J. R. Valverde ( 20 - aug - 1989):			*
*	    Version 2 terminated.				*
*	J. R. Valverde ( 08 - oct - 1989):			*
*	    Version 3 finished.					*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

/*  ----    General interest definitions    ----    */

#define	MAXLINE	    256		/* Longest input line			*/

#define	KEYS_NUMBER 1		/* Number of index keywords		*/

#define KEYLEN	    10		/* Length of keyword			*/

#define	PNAME_LEN   80		/* Maximum protein name length		*/

#define	MAX_SEQ_LEN 10000	/* Maximum sequence length		*/

#define MAX_REGIONS 100		/* Maximum number of regions 		*/

#define MAX_CLASSES 10		/* Maximum number of value categories	*/

#define AA_NUMBER   21		/* Number of aminoacids			*/

#define ASSIGNED    1		/* The region has already been used	*/

#define EXPLORED    1		/* This rule has already been checked	*/

#ifndef OK
#define	OK	    1		/* The protein passes this rule		*/
#endif

#define NOT_ALIGNED -1		/* Cannot align a region with a pattern */

#define F_INFINITY  9.99e+36	/* maximum float			*/

#define NEW(x)	    getmem( sizeof( x ) )/* safe space allocation macro	*/

#define P_EXT	    ".aa"	/* file extension for protein files	*/

#define D_EXT	    ".dxr"	/* file extension for descriptor files	*/

#define F_EXT	    ".nor"	/* default file extension for feat file	*/

/*  ----    Definitions for database access ----    */

#define DB_UNKNOWN	    0
#define DB_STANFORD	    1
#define DB_INTELLIGENETICS  2
#define DB_STRIDER	    3
#define DB_GENBANK	    4
#define DB_NBRF		    5
#define DB_PIR		    6
#define DB_PEARSON	    7
#define DB_EMBL		    8
#define DB_SWISSPROT	    9
#define DB_UWGCG	    10
#define DB_ZUKER	    11
#define DB_FITCH	    12
#define DB_PREER	    13

/*  ----    Definitions for queue manipulation	----	*/
                    
#define QUEUE_SIZE  128		/* size of static circular queue */

#define Q_FULL	    -1		/* error: queue is full	    */

#define Q_EMPTY	    NULL	/* error: queue is empty    */


/*  ---- Definitions for token manipulation while parsing files ----	*/

enum {
    NEW_RULE,
    DISTANCE,
    SEQUENCE,
    PROFILE,
    POSITION,
    RANGE,
    SEQ,
    IN_TOKEN,
    POS_RANGE,
    REGION,
    VALUES,
    VALUE,
    HYPHEN,
    EQUAL,
    LPAREN,
    RPAREN,
    COLON,
    POINTS,
    LESSTHAN
}   TOKENS;

#define NO_MORE_TOKENS	-1  /* No more tokens remain to be processed */


/*  ----    Definitions for program termination	    ----    */

/* Note: The following definitions should be modified according to the
system being used:
    EXIT_OK should be 0 for UNIX, 1 for VAX-VMS
    MAXPATH should be set to the appropriate value for the corresponding
	system.
*/

#ifdef VMS

/*
 *  Vax-VMS Condition codes (exception codes, not only exit codes):
 *
 *	31   28 27                16 15                 3 2   0
 *	+------+--------------------+--------------------+----+
 *	|      |                    |                    |    |
 *	+------+--------------------+--------------------+----+
 *	Control Facility number	     Message number	  Severity
 *
 *	Severity - Severity of the error condition. Bit 0 indicates
 *	success when set. Bits 1 & 2 distinguish degrees of success
 *	or failure:
 *	    0	Warning.    Execution continues. Unpredictable results
 *	    1	Success.    Execution continues. Expected results
 *	    2	Error.	    Execution continues. Erroneous results
 *	    3	Informat.   Execution continues. Information message
 *	    4	Fatal.	    Execution terminates. No output
 *	    5 - 7	    Reserved for DEC
 *	Message number - Number of the message associated with the
 *	error condition (the msg. may or not be displayed).
 *	Facility number - Number identifying the program originating
 *	the message. Bit 27 is set for user programs and clear for DEC's
 *	Control - Bit 28 inhibits the display of the error message. Bits
 *	29 - 31 are reserved for DEC.
 *
 * BIT_0    0x01
 * BIT_1    0x02
 * BIT_2    0x04
 * BIT_15   0x8000
 * BIT_16   0x10000
 * BIT_27   0x8000000
 * BIT_28   0x10000000
 *
 *		SevCode/User progr./Don't display
 */

#define VMS_WARNING (0 | 0x8000000 | 0x10000000)
#define VMS_SUCCESS (1 | 0x8000000 | 0x10000000)
#define VMS_ERROR   (2 | 0x8000000 | 0x10000000)
#define VMS_INFO    (3 | 0x8000000 | 0x10000000)
#define VMS_FATAL   (4 | 0x8000000 | 0x10000000)

#define EXIT_OK	    VMS_SUCCESS	    /* Successful termination	*/
#define BAD_ARGS    VMS_INFO	    /* Informative message	*/
#define LETAL_ERROR VMS_FATAL	    /* Letal error message	*/
#define	MAXPATH	    128		    /* Max length of file name	*/
#define OPTION_DELIMITER    '/'	    /* Options delimiter char	*/
#define DRIVECHAR	    ':'
#define BEGDIRCHAR	    '['	    /* use \ for MSDOS, : for Mac */
#define DIRDELIMITER	    '.'
#define ENDDIRCHAR	    ']'
#define FILEXTCHAR	    '.'

#else
#ifdef UNIX

#define UNIX_SUCCESS	0

#define EXIT_OK	    UNIX_SUCCESS
#define BAD_ARGS    1
#define LETAL_ERROR 2
#define MAXPATH	    128
#define OPTION_DELIMITER    '-'
#define DRIVECHAR	    ':'
#define BEGDIRCHAR	    '/'
#define DIRDELIMITER	    '/'
#define ENDDIRCHAR	    '/'
#define FILEXTCHAR	    '.'

#endif
#endif
                               
