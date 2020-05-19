# kill non-dlname lines
/^\(dlname\|libdir\)=/! { d }

/^dlname=/ {

# kill leading/trailing junk
s/^dlname='//

# kill from the last quote to the end
s/'.*$//

# kill blank lines
/./!d

# write out the lib on its own line
s/.*/\/&\n/g

# kill the EOL
s/\n$//

# hold it
h
}

/^libdir=/ {

# kill leading/trailing junk
s/^libdir='//

# kill from the last quote to the end
s/'.*$//

# paste
G

# kill the EOL
s/\n//
p
}