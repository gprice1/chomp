
import sys

def cleanList( list_object ):
    return filter( lambda x: len(x) > 3,
                   map( lambda x: x.strip(), list_object) )

def main():

    if len( sys.argv ) != 3:
        print "Usage: python foo.py <input_file> <output_file>"
        return

    fileReadName = sys.argv[1]
    fileWriteName = sys.argv[2]


    fileRead = open( fileReadName, "r" )
    fileWrite = open( fileWriteName, "w" )

    string = fileRead.read()

    lines = string.split( '\n' )
    lines = cleanList( lines )
    
    name_data = map( lambda x: x.split('{'), lines )
    name_data = map( lambda x: [x[0], x[1]], name_data )
                    
    print name_data[0][0], "-----", name_data[0][1]
    
    print name_data
    
    dictionary = {}
    
    for point in name_data:
        name = point[0]
        data = point[1]

        data_points = data.split( '[' )
        data_points = cleanList( data_points )
        data_points = map( lambda x: x.replace( ']', '' ), data_points )

        data_points = map( lambda x: x[:-1].replace(' ', '').split( "," ),
                           data_points)

        print name, "--", data_points[-1]

main()
