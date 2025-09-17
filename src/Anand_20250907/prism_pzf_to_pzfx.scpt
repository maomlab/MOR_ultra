on run argv
    if (count of argv) is not 2 then
        display dialog "Usage: osascript prism_pzf_to_pzfx.scpt /path/to/input.pzf /path/to/output.pzfx"
        return
    end if

    set inputFile to POSIX path of (item 1 of argv)
    set outputFile to POSIX path of (item 2 of argv)

    -- Validate extensions
    if inputFile does not end with ".pzf" then
        display dialog "Error: Input file must have .pzf extension"
        return
    end if
    if outputFile does not end with ".pzfx" then
        display dialog "Error: Output file must have .pzfx extension"
        return
    end if

    -- Generate sanitized temporary filenames
    set tmpBase to do shell script "mktemp -u /tmp/prismXXXXXX"
    set tmpInput to tmpBase & ".pzf"
    set tmpOutput to tmpBase & ".pzfx"
    set tmpScript to tmpBase & ".pzc"

    -- Ensure no leftover temp files
    do shell script "rm -f " & quoted form of tmpInput & " " & quoted form of tmpOutput & " " & quoted form of tmpScript

    -- Copy input to safe temporary file
    do shell script "cp " & quoted form of inputFile & " " & quoted form of tmpInput

    -- Build Prism script content
    set prismScript to "Open \"" & tmpInput & "\"\n" & ¬
                       "Save \"" & tmpOutput & "\"\n" & ¬
                       "Close\n"

    -- Write Prism script to temp file
    set tmpFile to open for access (POSIX file tmpScript as text) with write permission
    set eof of tmpFile to 0
    write prismScript to tmpFile
    close access tmpFile

    -- Run Prism with the script
    tell application "Prism 10"
        activate
        open POSIX file tmpScript
    end tell

    -- Wait a bit for Prism to finish (tune if needed for big projects)
    delay 3

    -- Verify that Prism produced the .pzfx
    set fileExists to false
    try
        do shell script "test -f " & quoted form of tmpOutput
        set fileExists to true
    end try

    if fileExists is false then
        display dialog "Error: Prism did not create the expected .pzfx file"
        -- Cleanup temp files
        do shell script "rm -f " & quoted form of tmpInput & " " & quoted form of tmpScript
        return
    end if

    -- Move the output back to requested location
    do shell script "mv -f " & quoted form of tmpOutput & " " & quoted form of outputFile

    -- Cleanup temp files
    do shell script "rm -f " & quoted form of tmpInput & " " & quoted form of tmpScript
end run
