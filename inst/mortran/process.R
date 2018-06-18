library(stringr)
library(tools)
##
## Jerry uses both plain fortran, escaped with %FORTRAN
## and mortran (mode forced with %MORTRAN after a %FORTRAN
## directive).
##
insertImplicitMortran <- function(inWhat = c("subroutine", "function"), f) {
    ## This handles the MORTRAN statements
    inWhat <- match.arg(inWhat)
    ## Find all definitions of inWhat at line beginning
    lineStart <- grep(paste0("^", inWhat),  f)
    lineEnd <- lineStart
    ## Since statements can span multiple lines
    ## figure out which is the last line of the
    ##  definition. I.e. look for semi-colon
    for (i in seq_along(lineEnd)) {
        x <- lineEnd[i];
        while(!grepl(";", str_trim(f[x])))  {
            x <- x + 1
        }
        lineEnd[i] <- x
    }
    ff <- as.list(f)
    ## Append implicit statement after each statement
    for (x in lineEnd) {
        ff[[x]] <- c(ff[[x]], "implicit double precision(a-h,o-z);")
    }
    unlist(ff)
}

insertImplicitFortran <- function(inWhat = c("subroutine", "function"), f) {
    ## This handles the FORTRAN statements
    inWhat <- match.arg(inWhat)
    padding = "      "  ## 6 spaces
    ## Find all definitions of inWhat at line beginning
    lineStart <- grep(paste0("^", padding, inWhat),  f)
    lineEnd <- lineStart
    ## Since statements can span multiple lines
    ## figure out which is the last line of the
    ##  definition. I.e. look for matched paren
    for (i in seq_along(lineEnd)) {
        x <- lineEnd[i];
        all_open_paren <- nrow(str_locate_all(x, "\\(")[[1]])
        all_close_paren <- nrow(str_locate_all(x, "\\)")[[1]])
        diff <- all_open_paren - all_close_paren
        while(diff != 0) {  ## assuming first open paren not continued!!
            x <- x + 1
            all_open_paren <- nrow(str_locate_all(x, "(")[[1]])
            all_close_paren <- nrow(str_locate_all(x, ")")[[1]])
            diff <- diff + (all_open_paren - all_close_paren)
        }
        lineEnd[i] <- x
    }
    ff <- as.list(f)
    ## Append implicit statement after each statement
    for (x in lineEnd) {
        ff[[x]] <- c(ff[[x]], ##paste0(padding, "implicit integer(i-n)"),
                     paste0(padding, "implicit double precision(a-h,o-z)"))
    }
    unlist(ff)
}

process_file <- function(inputMortranFile,
                         outputFile = paste0(file_path_sans_ext(basename(inputMortranFile)), ".f")) {

    outputMortranFile = paste0("Fixed-", file_path_sans_ext(basename(inputMortranFile)), ".m")
    cat("Processing Mortran; inserting implicit statements\n")
    f <- readLines(inputMortranFile)
    f <- insertImplicitMortran("subroutine", f) ## fix mortran subroutines
    f <- insertImplicitMortran("function", f)  ## fix mortran functions
    f <- insertImplicitFortran("subroutine", f) ## fix fortran subroutines
    f <- insertImplicitFortran("function", f)  ## fix fortran functions

    ## Next replace all real by double
    cat("Processing Mortran; replacing reals by double precision\n")
    f <- gsub("real", "double precision", f)

    ## Finally fix constants with e[+-]?[0-9]+.
    ##(?<first>[[:upper:]][[:lower:]]+) (?<last>[[:upper:]][[:lower:]]+)
    ##const.rex <- "([eE][+/]?[[:digit:]]+)"
    ##hits <- gregexpr(const.rex, f, perl= TRUE)

    cat("Checking for long lines; can cause problems downstream if not fixed\n")
    for (i in seq_along(f)) {
        x <- str_trim(f[i], "right")
        if (nchar(x) > 80) {
            cat(sprintf("%3d %s\n", i, x))
        }
    }

    writeLines(f, con = "temp.m")

    ## Now run mortran
    cat("Running Mortran\n")
    system2(command = "./m77", stdin = "temp.m")
    cat("Chopping Lines at 72 cols\n")
    code_lines <- substring(readLines(con = "./mo.for"), 1, 72)
    writeLines(code_lines, con = "temp.f")
    cat("Running gfortran to detect warning lines on unused labels\n")
    system2(command = "gfortran", args = c("-Wunused", "-c", "temp.f"), stderr = "gfortran.out")
    cat("Scanning gfortran output for warnings on unusued labels\n")
    warnings <- readLines("gfortran.out")
    line_numbers <- grep('^temp.f', warnings)
    label_warning_line_numbers <- grep(pattern = "^Warning: Label [0-9]+ at", warnings)
    ## Crude check that no other unused warnings besides labels are present
    nW <- length(label_warning_line_numbers)
    if (sum(label_warning_line_numbers[-nW] + 1 - line_numbers[-1]) > 0) {
        stop("Warnings besides numeric label warnings need fixing; stopping!")
    }
    for (i in seq_len(nW)) {
        offending_line <- as.integer(str_extract(warnings[line_numbers[i]], pattern="([0-9]+)"))
        code_line <- code_lines[offending_line]
        offending_label <- str_extract(warnings[label_warning_line_numbers[i]],
                                       pattern="([0-9]+)")
        code_lines[offending_line] <- sub(pattern = offending_label,
                                          replacement = str_pad("", width = nchar(offending_label)),
                                          x = code_lines[offending_line])
    }
    writeLines(code_lines, con = outputFile)
    file.rename(from = "temp.m", to = outputMortranFile)
    cat(sprintf("Fixed mortran in %s; fortran in %s\n", outputMortranFile, outputFile))
    invisible(TRUE)
}

process_file(inputMortranFile = "glasso.m")
