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
$MAIN = "main"; // Table name 
echo $username;
// Connect to server and select database. 
mysql_connect("$host", "$username", "$password")or die("cannot connect"); 
mysql_select_db($db_name)or die("cannot select DB: ". mysql_error()); 



// Get values from form 
$organism = mysql_escape_string($_POST['organism']);
$pathway = mysql_escape_string($_POST['pathway']); 
$linkout = mysql_escape_string($_POST['linkout']);
$dna_sequence = mysql_escape_string($_POST['dna_sequence']);
$uni_prot_id = mysql_escape_string($_POST['uni_prot_id']);
$norine_id = mysql_escape_string($_POST['norine_id']);

$invalid = false;
if (!$organism && !$pathway && !$linkout && !$dna_sequence && !$uni_prot_id && !$norine_id) {
    $clear_form = true;
} else {
    if (!$organism) {
        $organism_invalid = true;
        $invalid = true;
    }
    if (!$pathway) {
        $pathway_invalid = true;
        $invalid = true;
    }
    if (!$linkout) {
        $linkout_invalid = true;
        $invalid = true;
    }
    if (!$dna_sequence) {
        $dna_sequence_invalid = true;
        $invalid = true;
    }
    if (!$uni_prot_id) {
        $uni_prot_id_invalid = true;
        $invalid = true;
    }
    if (!$norine_id) {
        $norine_id_invalid = true;
        $invalid = true;
    }
}

if (!$invalid) {
             // Insert data into mysql 
$sql="INSERT INTO $MAIN SET
`organism` = '$organism',
`pathway` = '$pathway',
`linkout` = '$linkout',
`dna_sequence` = '$dna_sequence',
`uni_prot_id` = '$uni_prot_id',
`norine_id` ='$norine_id'";
$result = mysql_query($sql); 

// if successfully insert data into database, displays message "Successful". 
if($result){ 
echo "Successful"; 
echo "<BR>"; 
echo "<a href='insert.php'>Back to main page</a>";
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
        Organism: <input type="text" name="organism" class="<?php echo getClass($organism_invalid);?>" value="<?php echo getValue('organism');?>"/><br/>
        Pathway: <input type="text" name="pathway" class="<?php echo getClass($pathway_invalid);?>" value="<?php echo getValue('pathway');?>"/><br/>
        Linkout: <input type="text" name="linkout" class="<?php echo getClass($linkout_invalid);?>" value="<?php echo getValue('linkout');?>"/><br/>
        dna_sequence: <input type="text" name="dna_sequence" class="<?php echo getClass($dna_sequence_invalid);?>" value="<?php echo getValue('dna_sequence');?>"/><br/>
        UniProt ID: <input type="text" name="uni_prot_id" class="<?php echo getClass($uni_prot_id_invalid);?>" value="<?php echo getValue('uni_prot_id');?>"/><br/>
        NORINE ID: <input type="text" name="norine_id" class="<?php echo getClass($norine_id_invalid);?>" value="<?php echo getValue('norine_id');?>"/><br/>
        <input type="submit" value="Insert into database"/>
</form>


</body>
</html>
