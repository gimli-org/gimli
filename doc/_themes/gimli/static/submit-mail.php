<?php
if(isset($_POST['url']) && $_POST['url'] == ''){

    if(($_POST['name'] == '') |
        ($_POST['email'] == '') |
        ($_POST['subject'] == '') |
        ($_POST['message'] == '')) {
        header("Location: /doc/contact_failed.html");
        exit;
        }
    else {
	$to = 'mail@pygimli.org';
    $subject = $_POST['subject'];

    $name = $_POST['name'];
    $email = $_POST['email'];
    $body = "Message from ".$name." (sent via pygimli.org):\n\n";
    $body .= $_POST['message'];
    $from = "From: ".$email;
    mail($to, $subject, $body, $from);

    header("Location: /doc/contact_success.html");
    exit;
}}
?>
