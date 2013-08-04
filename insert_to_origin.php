<html>
<head>
    <title>Insert to table "origin"</title>
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
// open config-file
require("php_config.php");
$origin = "origin";

//connect to database
$connection = new mysqli("$host", "$username", "$password", $db_name);
if ($connection->connect_errno)
    die("cannot connect");

    // Get values from form
$organism = mysql_escape_string($_POST['organism']);
$pathway = mysql_escape_string($_POST['pathway']);
$linkout = mysql_escape_string($_POST['linkout']);
$dna_sequence = mysql_escape_string($_POST['dna_sequence']);
$uni_prot_id = mysql_escape_string($_POST['uni_prot_id']);
$norine_id = mysql_escape_string($_POST['norine_id']);
$description = mysql_escape_string($_POST['description']);

$invalid = false;
if (!$organism && !$pathway && !$dna_sequence && !$description) {
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
    if (!$dna_sequence) {
        $dna_sequence_invalid = true;
        $invalid = true;
    }
    if (!$description) {
        $description_invalid = true;
        $invalid = true;
    }
}

// Insert data into mysql
if (!$invalid) { 
$sql="INSERT INTO $origin SET
`organism` = '$organism',
`pathway` = '$pathway',
`linkout` = '$linkout',
`dna_sequence` = '$dna_sequence',
`uni_prot_id` = '$uni_prot_id',
`norine_id` = '$norine_id',
`description` ='$description'";
$result = mysql_query($sql);

// if successfully insert data into database, displays message "Successful". 
if($result){ 
echo "Successful"; 
echo "<BR>"; 
echo "<a href='insert4.php'>Back to main page</a>";
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
        Linkout: <input type="text" name="linkout" class="valid" value="<?php echo getValue('linkout');?>"/><br/>
        DNA Sequence: <input type="text" name="dna_sequence" class="<?php echo getClass($dna_sequence_invalid);?>" value="<?php echo getValue('dna_sequence');?>"/><br/>
        UniProt ID: <input type="text" name="uni_prot_id" class="valid" value="<?php echo getValue('uni_prot_id');?>"/><br/>
        Norine ID: <input type="text" name="norine_id" class="valid" value="<?php echo getValue('norine_id');?>"/><br/>
        Description: <input type="text" name="description" class="<?php echo getClass($description_invalid);?>" value="<?php echo getValue('description');?>" /><br/>
        <input type="submit" value="Insert into database"/>
</form>
</body>
</html>