/**
 * Created by Mada on 4/19/2018.
 */
$(document).ready(function () {
	
	$.ajax({
	  type: "POST",
	    url: "/cgi-bin/firstCgi.py",
	      data: {"test":"test"},
	      dataType: "json",
	        success: function(result){
			console.log(result);
		},
		
		  });
});
