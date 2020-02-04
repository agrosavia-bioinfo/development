# JAVA_LD_LIBRARY_PATH library path necessary at run-time
# JAVA_CPPFLAGS  C preprocessor flags necessary to compile JNI programs
# JAVA_LIBS      libraries (as linker flags) necessary to compile JNI programs

export JAVA_HOME=/opt/apps/jdk
export JAVA_LD_LIBRARY_PATH=$JAVA_HOME/jre/lib/amd64/server
export JAVA_CPPFLAGS="-I$JAVA_HOME/include -I$JAVA_HOME/include/linux "
export JAVA_LIBS="-L$JAVA_HOME/jre/lib/amd64/server -ljvm" 

env|grep JAVA

R CMD javareconf

