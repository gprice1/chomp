
import sys

def cleanList( list_object ):
    return filter( lambda x: len(x) > 3,
                   map( lambda x: x.strip(), list_object) )

def parseToString( name, data, list_of_fields ):
    
    string = name + "\n"

    for field in list_of_fields:
        
        for point in data:
            string += point[field] + "\t"
            
        string += "\n"
    return string

def parseAllData( data, list_of_fields ):
    big_string = ""
    
    for point in data:
        name = point[0]
        info = point[1]
        
        string = parseToString( name, info[1:], list_of_fields )

        big_string += string

    return big_string


def main():

    if len( sys.argv ) != 3:
        print "Usage: python foo.py <input_file> <output_file>"
        return

    fileReadName = sys.argv[1]
    fileWriteName = sys.argv[2]


    fileRead = open( fileReadName, "r" )
    fileWrite = open( fileWriteName, "w" )

    string = fileRead.read()
    string = string.replace( " ", "" )
    string = string.replace( "}", "" )

    lines = string.split( '\n' )
    lines = cleanList( lines )
    
    name_data = map( lambda x: x.split('{'), lines )
    name_data = map( lambda x: [x[0], x[1]], name_data )
                    
    
    string = ""
    for point in name_data:
        name = point[0]
        data = point[1]

        data = data.split( '[' )
        data = cleanList( data )
        data = map( lambda x: x.replace( ']', '' ), data )
        data = map( lambda x: x.split( "," ), data)
        data = filter( lambda x: len(x) > 0, data )
        data = [ filter( lambda x: len(x) > 0, y ) for y in data ]
        data = [ map( lambda x: x.split(':'), y ) for y in data ]
        data = [ dict( y ) for y in data ]

        
        print name, data[-1]
        string += parseToString( name, data[1:], ['total', 'objective'] )

    fileWrite.write( string )

main()
