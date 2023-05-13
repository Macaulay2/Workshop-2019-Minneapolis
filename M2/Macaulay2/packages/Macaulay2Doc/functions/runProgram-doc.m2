doc ///
  Key
    ProgramRun
  Headline
    result of running an external program
  Description
    Text
      A hash table returned by @TO runProgram@ with the following
      strings as keys:

      @UL {
	 {TT "\"command\"", ", the command that was used to run the program."},
	 {TT "\"output\"", ", the output of the program to ", TT "stdout", "."},
	 {TT "\"error\"", ", the output of the program to ", TT "stderr", "."},
	 {TT "\"return value\"", ", the return value of the program, ",
	     "possibly multiplied by 256 (see ", TO "run", ").  ",
	     "Note that this is what is displayed when printing a ",
	     TT "ProgramRun", " object."}}@

      In addition, if @TO runProgram@ is called with the @TT "KeepFiles"@
      option set to @TO true@, then the following keys will be present:

      @UL {
	 {TT "\"output file\"", ", the path to a file containing the output ",
	     "of the program to ", TT "stdout", "."},
	 {TT "\"error file\"", ", the path to a file containing the output ",
	     "of the program to ", TT "stderr", "."}}@
  SeeAlso
    runProgram
///

doc ///
  Key
    KeepFiles
  Headline
    whether to keep temporary files
///

doc ///
  Key
    RunDirectory
  Headline
    the directory from which to run a program
///

doc ///
  Key
    runProgram
    (runProgram, Program, String)
    (runProgram, Program, String, String)
    (runProgram, String, String)
    [runProgram, KeepFiles]
    [runProgram, RaiseError]
    [runProgram, RunDirectory]
    [runProgram, Verbose]
  Headline
    run an external program
  Usage
    runProgram(program, args)
    runProgram(program, exe, args)
  Inputs
    program:Program
      the program to run, generated by @TO findProgram@.
    exe:String
      the specific executable file to run. This is only necessary if
      the program consists of multiple such files.  If not given, then
      @TT "program#\"name\""@ is used.
    args:String
      the command line arguments passed to the program.
    KeepFiles => Boolean
      whether to keep the temporary files containing the program's output.
    RaiseError => Boolean
      whether to raise an error if the program returns a nonzero value.
    RunDirectory => String
      the directory from which to run the program.  If it does not
      exist, then it will be created.  If @TO null@, then the program
      will be run from the current working directory (see
      @TO currentDirectory@).
    Verbose => Boolean
      whether to print the command line input and the program's output.
  Description
    Text
      This method runs an external program which has already been
      loaded using @TO findProgram@.  The results of this run are
      available in a @TO ProgramRun@ object.
    Example
      gfan = findProgram("gfan", "gfan --help")
      runProgram(gfan, "_version")
      oo#"output"
      runProgram(gfan, "_foo", RaiseError => false)
      oo#"error"
    Text
      The value corresponding to the @TT "\"output\""@ key may also be
      obtained using @TO toString@.
    Example
      toString runProgram(gfan, "_version")
    Text
      It is also possible to skip @TO findProgram@ and just provide
      two strings: the name of the program and the command line
      arguments.
    Example
      runProgram("normaliz", "--version")
      peek oo
  SeeAlso
    ProgramRun
    findProgram
    (status, ProgramRun)
///

doc ///
  Key
    (status, ProgramRun)
  Headline
    get the return status of a program run
  Usage
    status pr
  Inputs
    pr:ProgramRun
  Outputs
    :ZZ
  Description
    Text
      Get the return status of a program run.  Usually, 0 means that it was
      successful.
    Example
      normaliz = findProgram "normaliz"
      status runProgram(normaliz, "--version")
  SeeAlso
    findProgram
    runProgram
///

doc ///
  Key
    (symbol <<, Program, Thing)
  Headline
    run program with input redirection
  Usage
    prog << x
  Inputs
    prog:Program
    x:Thing
  Outputs
    :ProgramRun
  Description
    Text
      Write @TT "x"@ to a temporary file and run @TT "prog"@ with this
      file as input using input redirection (the @TT "<"@ operator in
      a POSIX shell).
    Example
      M2 = findProgram "M2"
      M2 << "2 + 2"
      toString oo
  SeeAlso
    findProgram
    runProgram
///