<html>
<head>
    <title>Insert NRPS</title>
    <style type="text/css">
        .valid {
            background-color: #FFFFFF;
        }
        .invalid {
            background-color: #FF0000;
        }
    </style>
</head>

<body>
    
<?php 

$host = "localhost"; // Host name 
$username = "root"; // Mysql username 
$password = ""; // Mysql password 
$db_name = "nrps_designer"; // Database name 
$substrates = "substrates"; // Table name

echo $username;

// Connect to server and select database. 
mysql_connect("$host", "$username", "$password")or die("cannot connect"); 
mysql_select_db($db_name)or die("cannot select DB: ". mysql_error());

//values in form
$aminoacid = mysql_escape_string($_POST['aminoacid']);
$linkout = mysql_escape_string($_POST['linkout']); 
$structure = mysql_escape_string($_POST['structure']);

$invalid = false;
if (!$aminoacid && !$linkout && !$structure ) {
    $clear_form = true;
} else {
    if (!$aminoacid) {
        $aminoacid_invalid = true;
        $invalid = true;
    }
    if (!$linkout) {
        $linkout_invalid = true;
        $invalid = true;
    }
    if (!$structure) {
        $structure_invalid = true;
        $invalid = true;
    }
}

if (!$invalid) {   
             // Insert data into mysql 
$sql="INSERT INTO $substrates SET
`aminoacid` = '$aminoacid',
`linkout` = '$linkout',
`structure` = '$structure'";
$result = mysql_query($sql);

//if successfull insert data into database
if($result){ 
echo "Successful"; 
echo "<BR>"; 
echo "<a href='insert2.php'>Back to main page</a>";
$clear_form = true;
} else { 
    echo mysql_error(); 
} 
    
}

// close connection 
mysql_close();
function getClass($invalid) {
    return ($invalid ? "invalid" : "valid");
}
function getValue($var) {
    global $clear_form;
    return ($clear_form ? "" : $_POST[$var]);
}
?>
<form method="post" action="<?php echo $_SERVER["PHP_SELF"];?>">
        Aminoacid: <input type="text" name="aminoacid" class="<?php echo getClass($aminoacid_invalid);?>" value="<?php echo getValue('aminoacid');?>"/><br/>
        Linkout: <input type="text" name="linkout" class="<?php echo getClass($linkout_invalid);?>" value="<?php echo getValue('linkout');?>"/><br/>
        Structure: <input type="text" name="structure" class="<?php echo getClass($structure_invalid);?>" value="<?php echo getValue('structure');?>"/><br/>
        <input type="submit" value="Insert into database"/>
</form>


</body>
</html>
