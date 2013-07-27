-- phpMyAdmin SQL Dump
-- version 3.4.11.1deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Jul 16, 2013 at 05:27 PM
-- Server version: 5.5.31
-- PHP Version: 5.4.6-1ubuntu1.2

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";

--
-- Database: `nrps_designer`
--

CREATE DATABASE `nrps_designer` DEFAULT CHARACTER SET utf8 COLLATE utf8_unicode_ci;
USE `nrps_designer`;

-- --------------------------------------------------------
--
-- Table structure for table `domains`
--

DROP TABLE IF EXISTS `domains`;
CREATE TABLE IF NOT EXISTS `domains` (
  `domain_id` int(20) NOT NULL,
  `module_id` int(20) NOT NULL,
  `pathway_id` int(20) NOT NULL,
  `refseq_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `type` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `substrate_specificity_aa_id` int(20),
  `chirality` tinyint(1) NOT NULL,
  `uni_prot_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`domain_id`),
  KEY `substrate_specificity_aa_id` (`substrate_specificity_aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Table structure for table `main`
--

DROP TABLE IF EXISTS `main`;
CREATE TABLE IF NOT EXISTS `main` (
  `pathway_id` int(20) NOT NULL AUTO_INCREMENT,
  `organism` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `pathway` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `linkout` varchar(300) COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `uni_prot_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `norine_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=4 ;

-- --------------------------------------------------------

--
-- Table structure for table `substrates`
--

DROP TABLE IF EXISTS `substrates`;
CREATE TABLE IF NOT EXISTS `substrates` (
  `aminoacid` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `aa_id` int(20) NOT NULL AUTO_INCREMENT,
  `structure` text COLLATE utf8_unicode_ci NOT NULL,
  `modification` int(20),
  `linkout` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Constraints for dumped tables
--

--
-- Table structure for table `modifications`
--

DROP TABLE IF EXISTS `modifications`;
CREATE TABLE IF NOT EXISTS `modifications` (
  `type` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `mod_id` int(20) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`mod_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;



--
-- Constraints for tables `domain_types`
-- reference to main table
ALTER TABLE `domains`
  ADD CONSTRAINT `domains_ibfk_1` FOREIGN KEY (`pathway_id`) REFERENCES `main` (`pathway_id`) ON DELETE CASCADE ON UPDATE CASCADE;
  
ALTER TABLE `domains`
  ADD INDEX ( `substrate_specificity_aa_id` ),
  ADD FOREIGN KEY (`substrate_specificity_aa_id`) REFERENCES `substrates` (`aa_id`) ON DELETE SET NULL ON UPDATE CASCADE;

