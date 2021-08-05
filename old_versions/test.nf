//String pattern = "MMddyy";
//SimpleDateFormat simpleDateFormat = new SimpleDateFormat(pattern);
x = new java.util.Date().format( 'MMddyy' )
println x

process get_changes {
    output:
    path("lazy_map.txt") into map_ch
    
    when:
    ( params.lazy )

    shell:
    
    '''
    echo "yes" > "lazy_map.txt"
    '''
}